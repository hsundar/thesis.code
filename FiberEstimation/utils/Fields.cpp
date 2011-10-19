

/**
 * \mainpage Contains code for various operations such as smoothing, etc. on tensor, vector and scalar fields
 * Supports several executables for processing diffusion tensor data.
 */

# include "Fields.h"
#include <iostream>
using namespace std;



//------
void ScalarField::ComputeSmoothingFilter(float sigma)
//------
{  
  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  float total_sm= 0.0;

  if(Xdim%2 ==0) { std::cout<<"ComputeSmoothingFilter: Kernel size "<<Xdim<<" should be odd\n"; exit(0); }
  if(Ydim%2 ==0) { std::cout<<"ComputeSmoothingFilter: Kernel size "<<Ydim<<" should be odd\n"; exit(0); }
  if(Zdim%2 ==0) { std::cout<<"ComputeSmoothingFilter: Kernel size "<<Zdim<<" should be odd\n"; exit(0); }

  // assume all sizes are odd
  int cen_x = (Xdim-1)/2;
  int cen_y = (Ydim-1)/2;
  int cen_z = (Zdim-1)/2; 
 
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	 
	  float x= (ix-cen_x)*Xres;
	  float y= (iy-cen_y)*Yres;
	  float z= (iz-cen_z)*Zres;
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  m_Data[ind] = exp(-(x*x + y*y + z*z)/(2*sigma*sigma));

	  total_sm += m_Data[ind];

	}
    
  // normalize
  
  for(iz=0; iz<Zdim; iz++)
    for(ix=0; ix<Xdim; ix++)
      {
	for(iy=0; iy<Ydim; iy++)
	  {
	    int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	    m_Data[ind] = m_Data[ind]/total_sm;
	    //  if(iz==Zdim/2) std::cout<<m_Data[ind]<< " ";
	  }
	//	std::cout<<"\n";
      }
  //std::cout<<"Smoothing filter of width "<<sigma<<" computed\n";

}

//------
void ScalarField::ComputeSmoothingFilter(float sigma, float cen_x, float cen_y, float cen_z)
//------
{  
  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  float total_sm= 0.0;

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	 
	  float x= (ix-cen_x)*Xres;
	  float y= (iy-cen_y)*Yres;
	  float z= (iz-cen_z)*Zres;
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  m_Data[ind] = exp(-(x*x + y*y + z*z)/(2*sigma*sigma));

	  total_sm += m_Data[ind];

	}
    
  // normalize
  
  for(iz=0; iz<Zdim; iz++)
    for(ix=0; ix<Xdim; ix++)
      {
	for(iy=0; iy<Ydim; iy++)
	  {
	    int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	    m_Data[ind] = m_Data[ind]/total_sm;
	    //  if(iz==Zdim/2) std::cout<<m_Data[ind]<< " ";
	  }
	//	std::cout<<"\n";
      }
  //std::cout<<"Smoothing filter of width "<<sigma<<" computed\n";

}



//----
void ScalarField::Smooth(ScalarField *smoothingFilter, ScalarField *out)
//----
{

  // should consider doing this computation in the Fourier domain

  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  

  int ix1, iy1, iz1;
  int ix2, iy2, iz2;

  int fXdim = smoothingFilter->getSize(0);
  int fYdim = smoothingFilter->getSize(1);
  int fZdim = smoothingFilter->getSize(2);

  /*--
  TensorField dti_out;
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  --*/

  if(fXdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }
  if(fYdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }
  if(fZdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }

  // assume kernel sizes are odd
  
  float *kernel_Data = smoothingFilter->getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  for(ix=0; ix<Xdim; ix++)
    {
      for(iy=0; iy<Ydim; iy++)
	for(iz=0; iz<Zdim; iz++)
	  {
	  
	 
	  float *val= out->getAt(ix,iy,iz);
	  val[0] = 0.0;

	  if(ix < order_x || ix >= Xdim- order_x || iy<order_y || iy >= Ydim- order_y || iz<order_z || iz>=Zdim- order_z )
	    continue;

	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;

		  ix2 = ix -ix1 + order_x;
		  iy2 = iy -iy1 + order_y;
		  iz2 = iz -iz1 + order_z;
		  
		  int ind = iz2*Xdim*Ydim + iy2*Xdim + ix2;
		  		
		  val[0] += kernel_Data[ind1]*m_Data[ind];

     		}
	  }
      //std::cout<<ix<<" YZ slice smoothing done\n";
    }
  
}


//-------------------------------------------------------------

//----
void VectorField::AddOffset()
//----
{
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix,iy,iz,count;
  float *vec_data = getAt(0,0,0);
  std::cout<<"Adding offset\n";
  ix=Xdim/2; iy=Ydim/2; iz= Zdim/2;
  count= iz*Xdim*Ydim + iy*Xdim + ix;
  std::cout<<"\nCentral Original Vector: "<<vec_data[count*3]<<" "<<vec_data[count*3+1]<<" "<<vec_data[count*3+2]<<"\n";
  count=0; 
  
  for(iz=0; iz<Zdim; iz++)
    for(iy=0; iy<Ydim; iy++)
      for(ix=0; ix<Xdim; ix++)
	{
	  vec_data[count*3+0] = vec_data[count*3+0] + ix;
	  vec_data[count*3+1] = vec_data[count*3+1] + iy;
	  vec_data[count*3+2] = vec_data[count*3+2] + iz;
	  count++;
	} 
  ix=Xdim/2; iy=Ydim/2; iz= Zdim/2;
  count= iz*Xdim*Ydim + iy*Xdim + ix;
  std::cout<<"Central Modified Vector: "<<vec_data[count*3]<<" "<<vec_data[count*3+1]<<" "<<vec_data[count*3+2]<<"\n";
  
}

//----
void VectorField::SubtractOffset()
//----
{
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix,iy,iz,count;
  float *vec_data = getAt(0,0,0);
  std::cout<<"Subtracting offset\n";
  ix=Xdim/2; iy=Ydim/2; iz= Zdim/2;
  count= iz*Xdim*Ydim + iy*Xdim + ix;
  std::cout<<"\nCentral Original Vector: "<<vec_data[count*3]<<" "<<vec_data[count*3+1]<<" "<<vec_data[count*3+2]<<"\n";
  count=0; 
  
  for(iz=0; iz<Zdim; iz++)
    for(iy=0; iy<Ydim; iy++)
      for(ix=0; ix<Xdim; ix++)
	{
	  vec_data[count*3+0] = vec_data[count*3+0] - ix;
	  vec_data[count*3+1] = vec_data[count*3+1] - iy;
	  vec_data[count*3+2] = vec_data[count*3+2] - iz;
	  count++;
	} 
  ix=Xdim/2; iy=Ydim/2; iz= Zdim/2;
  count= iz*Xdim*Ydim + iy*Xdim + ix;
  std::cout<<"Central Modified Vector: "<<vec_data[count*3]<<" "<<vec_data[count*3+1]<<" "<<vec_data[count*3+2]<<"\n";
  
}

//----
void VectorField::FlipComp(int FlipDirn)
//----
{
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix,iy,iz,count, j;
  count=0;
  for(iz=0; iz<Zdim; iz++)
    for(iy=0; iy<Ydim; iy++)
      for(ix=0; ix<Xdim; ix++)
	{
	  m_Data[count*3+FlipDirn] = - m_Data[count*3+FlipDirn];
	  count++;
	}
}

//----
void VectorField::Normalize()
//----
{
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix,iy,iz,count,j;
  count=0;
  for(iz=0; iz<Zdim; iz++)
    for(iy=0; iy<Ydim; iy++)
      for(ix=0; ix<Xdim; ix++)
	{
	  float VecNorm= 0.0;

	  for(j=0; j<3; j++)
	    VecNorm += m_Data[count*3+j]*m_Data[count*3+j];

	  VecNorm = sqrt(VecNorm);

	  for(j=0; j<3; j++)
	    m_Data[count*3+j] /= VecNorm;

	  count++;
	}
}

// takes a forward mapping and returns a reverse mapping
//-----
void VectorField::reverseWarpField(VectorField *out)
//-----
{ 
  // starts 
  // estimate reverse vector field
  // this is a simple reverser : a more complex one is implemented by DGS
  // unfortunately, this simple one does not follow HAMMER coordinate conventions
  // reverseWarpFieldHAMMER follows these conventions

  float eps = 1e-5;

  // assume fwd field is from subject to template

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  int ix,iy,iz,count,j;
  int ix_template,iy_template,iz_template;

  ScalarField FwdVoxelCounts;
  
  FwdVoxelCounts.init(Xdim,Ydim,Zdim);
  FwdVoxelCounts.setVoxelSize(0,Xres); 
  FwdVoxelCounts.setVoxelSize(1,Yres); 
  FwdVoxelCounts.setVoxelSize(2,Zres);
  
  // initialize to zero
  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  
	  float *val= out->getAt(ix,iy,iz);	  
	  for(j=0; j < 3; j++)
	    val[j] = 0.0f;

	  val= FwdVoxelCounts.getAt(ix,iy,iz);
	  val[0] = 0.0f;
	  
	}

  int mutliply_mapped_voxels= 0;
  int unmapped_voxels= 0;
  int outside_voxels = 0;
  int singly_mapped_voxels= 0;


  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{

	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *warpField_val = getAt(ix,iy,iz);
	  ix_template = rint(warpField_val[0]);
	  iy_template = rint(warpField_val[1]);
	  iz_template = rint(warpField_val[2]);
	  
	  if(ix_template < 0 || iy_template < 0 || iz_template < 0 || ix_template >= Xdim || iy_template >= Ydim || iz_template >= Zdim )
	    {	     
	      outside_voxels++;
	      continue;	      
	    }
	  
	  float *rev_warpField_val = out->getAt(ix_template,iy_template,iz_template);
	  rev_warpField_val[0] += ix;
	  rev_warpField_val[1] += iy;
	  rev_warpField_val[2] += iz;
	  
	  float *FwdVoxelCounts_val = FwdVoxelCounts.getAt(ix_template,iy_template,iz_template);
	  FwdVoxelCounts_val[0] += 1.0f;
	  
	}

  std::cout<<"Step 1 (mapping) done "<<"\n";
  std::cout<<"Outside voxels count: "<<outside_voxels<<"\n";

  // average the displacements when >1 voxel in subject gets mapped to one voxel in template
  
   for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *FwdVoxelCounts_val = FwdVoxelCounts.getAt(ix,iy,iz);
	  if(FwdVoxelCounts_val[0]>1+eps) //>1
	    {
	      mutliply_mapped_voxels++;
	      float *rev_warpField_val = out->getAt(ix,iy,iz);
	      for(int tempi=0; tempi<3; tempi++)
		rev_warpField_val[tempi] /= FwdVoxelCounts_val[0];
	      continue;
	    }
	  
	  if(FwdVoxelCounts_val[0]>eps) //==1
	    {
	      singly_mapped_voxels++;
	      continue;
	    }
	  
	  if(FwdVoxelCounts_val[0]>-1+eps) //==0
	    {
	      unmapped_voxels++;
	      float *rev_warpField_val = out->getAt(ix,iy,iz);
	      rev_warpField_val[0] = ix;
	      rev_warpField_val[1] = iy;
	      rev_warpField_val[2] = iz;
	      continue;
	    }
	  
	}

   std::cout<<"Step 2 (averaging) done\n";
   
   std::cout<<"\nMultiply mapped voxels count: "<<mutliply_mapped_voxels<<"\nSingly mapped voxels count: "<<singly_mapped_voxels<<"\nUnmapped voxels count: "<<unmapped_voxels<<"\n";

   if(mutliply_mapped_voxels+ singly_mapped_voxels+ unmapped_voxels != Xdim*Ydim*Zdim)
     std::cout<<"\nError in voxel processing in computing reverse vector field\n";

}  

// takes a forward displacement field and returns a reverse displacement field
//-----
void VectorField::reverseWarpFieldHAMMER(VectorField *out)
//-----
{ 
  // estimate reverse vector field
  // this is a simple reverser : a more complex one is implemented by DGS

  float eps = 1e-5;

  // assume fwd field is from subject to template

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  int ix,iy,iz,count,j;
  int ix_template,iy_template,iz_template;

  ScalarField FwdVoxelCounts;
  
  FwdVoxelCounts.init(Xdim,Ydim,Zdim);
  FwdVoxelCounts.setVoxelSize(0,Xres); 
  FwdVoxelCounts.setVoxelSize(1,Yres); 
  FwdVoxelCounts.setVoxelSize(2,Zres);
  
  // initialize to zero
  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  
	  float *val= out->getAt(ix,iy,iz);	  
	  for(j=0; j < 3; j++)
	    val[j] = 0.0f;

	  val= FwdVoxelCounts.getAt(ix,iy,iz);
	  val[0] = 0.0f;
	  
	}

  int mutliply_mapped_voxels= 0;
  int unmapped_voxels= 0;
  int outside_voxels = 0;
  int singly_mapped_voxels= 0;

  // field from subject 2 template

  // (ix,iy,iz) now in subject

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  
	  float *warpField_val = getAt(ix,iy,iz);

	  // note that the (1,0) switch is needed to deal with HAMMER conventions

	  ix_template = rint(warpField_val[1] + ix);
	  iy_template = rint(warpField_val[0] + iy);
	  iz_template = rint(warpField_val[2] + iz);
	  
	  if(ix_template < 0 || iy_template < 0 || iz_template < 0 || ix_template >= Xdim || iy_template >= Ydim || iz_template >= Zdim )
	    {	     
	      outside_voxels++;
	      continue;	      
	    }
	  
	  float *rev_warpField_val = out->getAt(ix_template,iy_template,iz_template);
	  rev_warpField_val[0] += ix;
	  rev_warpField_val[1] += iy;
	  rev_warpField_val[2] += iz;
	  
	  float *FwdVoxelCounts_val = FwdVoxelCounts.getAt(ix_template,iy_template,iz_template);
	  FwdVoxelCounts_val[0] += 1.0f;
	  
	}

  std::cout<<"Step 1 (mapping) done "<<"\n";
  std::cout<<"Outside voxels count: "<<outside_voxels<<"\n";

  // average the displacements when >1 voxel in subject gets mapped to one voxel in template
  
   for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *FwdVoxelCounts_val = FwdVoxelCounts.getAt(ix,iy,iz);
	  if(FwdVoxelCounts_val[0]>1+eps) //>1
	    {
	      mutliply_mapped_voxels++;
	      float *rev_warpField_val = out->getAt(ix,iy,iz);
	      for(int tempi=0; tempi<3; tempi++)
		rev_warpField_val[tempi] /= FwdVoxelCounts_val[0];
	      continue;
	    }
	  
	  if(FwdVoxelCounts_val[0]>eps) //==1
	    {
	      singly_mapped_voxels++;
	      continue;
	    }
	  
	  if(FwdVoxelCounts_val[0]>-1+eps) //==0
	    {
	      unmapped_voxels++;
	      float *rev_warpField_val = out->getAt(ix,iy,iz);
	      rev_warpField_val[0] = ix;
	      rev_warpField_val[1] = iy;
	      rev_warpField_val[2] = iz;
	      continue;
	    }
	  
	}

   std::cout<<"Step 2 (averaging) done\n";

   // return to hammer convention

   for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *rev_warpField_val = out->getAt(ix,iy,iz);
	  float tmpx= rev_warpField_val[0];
	  float tmpy= rev_warpField_val[1];
	  float tmpz= rev_warpField_val[2];
	  rev_warpField_val[0]= tmpy - iy;
	  rev_warpField_val[1]= tmpx - ix;
	  rev_warpField_val[2]= tmpz - iz;
	}

   std::cout<<"\nMultiply mapped voxels count: "<<mutliply_mapped_voxels<<"\nSingly mapped voxels count: "<<singly_mapped_voxels<<"\nUnmapped voxels count: "<<unmapped_voxels<<"\n";

   if(mutliply_mapped_voxels+ singly_mapped_voxels+ unmapped_voxels != Xdim*Ydim*Zdim)
     std::cout<<"\nError in voxel processing in computing reverse vector field\n";

}  

//----
void VectorField::Smooth(ScalarField *smoothingFilter, VectorField *out)
//----
{

  // should consider doing this computation in the Fourier domain

  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  

  int ix1, iy1, iz1;
  int ix2, iy2, iz2;

  int fXdim = smoothingFilter->getSize(0);
  int fYdim = smoothingFilter->getSize(1);
  int fZdim = smoothingFilter->getSize(2);

  /*--
  TensorField dti_out;
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  --*/

  if(fXdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }
  if(fYdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }
  if(fZdim%2 ==0) { std::cout<<"Smooth: Kernel size should be odd\n"; exit(0); }

  // assume kernel sizes are odd
  
  float *kernel_Data = smoothingFilter->getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  for(ix=0; ix<Xdim; ix++)
    {
      for(iy=0; iy<Ydim; iy++)
	for(iz=0; iz<Zdim; iz++)
	  {
	  
	 
	  float *val= out->getAt(ix,iy,iz);
	  for(j=0; j<3; j++)
	    val[j] = 0.0;

	  if(ix < order_x || ix >= Xdim- order_x || iy<order_y || iy >= Ydim- order_y || iz<order_z || iz>=Zdim- order_z )
	    continue;

	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;

		  ix2 = ix -ix1 + order_x;
		  iy2 = iy -iy1 + order_y;
		  iz2 = iz -iz1 + order_z;
		  
		  int ind = iz2*Xdim*Ydim + iy2*Xdim + ix2;

		  for(j=0; j<3; j++)
		    val[j] += kernel_Data[ind1]*m_Data[ind*3 + j];

     		}
	  }
      //std::cout<<ix<<" YZ slice smoothing done\n";
    }
  
}


//-----
int VectorField::getVectorAtAnyPosition(float px, float py, float pz, struct Pt3d & P)
//-----
//---- Function directly taken from XDR's code

{
  // these 4 lines are added
  int dimZ= m_Dims[2];
  int dimX= m_Dims[0];
  int dimY= m_Dims[1];
  int planSize = dimX*dimY;
  
//Tri-linear interpolation in *fdMap according to the coordinate of smpt
 //This volume have sampZ slices in Z direction
 //This function is revised from  class VolumeData::getValue(fPoint pt)
 //  in VData.C
 // return 1 when P is properly assigned
 // return 0 when P is undefined (but set to be a 0-vector)
    float vx,vy,vz;
    int row,col,slice;
    struct Pt3d *pp, *pp1,*tpp,*tpp1;
    struct Pt3d v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;
    //struct Pt3d* volume=fdMap;
    struct Pt3d* volume=(struct Pt3d* ) m_Data; // only line modified
    int sampZ = dimZ;

    vx=px;vy=py;vz=pz;

	//cout <<"   BBB:  AnyPos: ("<<vx<<","<<vy<<","<<vz<<")"<<endl;
    if ( ((vx<0)||(vx>dimX-1)) || ((vy<0)||(vy>dimY-1)) 
        || ((vz<0)||(vz>sampZ-1)))
       return  0; //means value undefined

    col = int(vx); //cout << "col=" << col;
    row = int(vy); //cout << ",row=" << row;
    slice=int(vz); //cout << ",slice=" << slice << endl;
   
    pp = volume + planSize*slice + row*dimX + col;
    if (slice == sampZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = *((struct Pt3d*)pp); //cout << "v1=" << v1.x <<","<< v1.y<<"," <<v1.z <<endl;
    v5 = *((struct Pt3d*)pp1);//cout << "v5=" << v5.x<<","<<v5.y<<","<<v5.z << endl;

    if (col == dimX-1){
	v2 = v1;
	v6 = v5;
    } else {
        v2 = *((struct Pt3d*)(pp+1));
    	v6 = *((struct Pt3d*)(pp1+1));
    }
    if (row == dimY-1) {
	v3 = v1;tpp=pp;
	v7 = v5;tpp1=pp1;
    } else {
    	tpp=pp+dimX;   v3 = *((struct Pt3d*)(tpp));
    	tpp1=pp1+dimX; v7 = *((struct Pt3d*)(tpp1));
    }
    if (col == dimX-1){
	v4 = v3;
    	v8 = v7;
    } else {
	v4 = *((struct Pt3d*)(tpp+1));
    	v8 = *((struct Pt3d*)(tpp1+1));
    } 
    //tri-linear interpolation
    tx = vx - col; ty = vy - row; tz = vz - slice;
        /*
    	cout << "v2=" << v2.x << ","<< v2.y <<","<<v2.z <<endl;
      	cout << "v6=" << v6.x << ","<< v6.y <<","<<v6.z <<endl;
	cout << "v4=" << v4.x << ","<< v4.y <<","<<v4.z <<endl;
	cout << "v8=" << v8.x << ","<< v8.y <<","<<v8.z <<endl;

	cout << "tx=" << tx << endl;
	cout << "ty=" << ty << endl;
	cout << "tz=" << tz << endl;*/
	
    //v1 & v3
    v1.x = tx*v2.x+(1-tx)*v1.x;
    v1.y = tx*v2.y+(1-tx)*v1.y;
    v1.z = tx*v2.z+(1-tx)*v1.z; 
    v3.x = tx*v4.x+(1-tx)*v3.x;
    v3.y = tx*v4.y+(1-tx)*v3.y;
    v3.z = tx*v4.z+(1-tx)*v3.z;

    // v5 & v7
    v5.x = tx*v6.x+(1-tx)*v5.x;
    v5.y = tx*v6.y+(1-tx)*v5.y;
    v5.z = tx*v6.z+(1-tx)*v5.z;
    v7.x = tx*v8.x+(1-tx)*v7.x;
    v7.y = tx*v8.y+(1-tx)*v7.y;
    v7.z = tx*v8.z+(1-tx)*v7.z;

    v1.x = ty*v3.x+(1-ty)*v1.x;
    v1.y = ty*v3.y+(1-ty)*v1.y;
    v1.z = ty*v3.z+(1-ty)*v1.z;
    v5.x = ty*v7.x+(1-ty)*v5.x;
    v5.y = ty*v7.y+(1-ty)*v5.y;
    v5.z = ty*v7.z+(1-ty)*v5.z;

    v1.x = tz*v5.x+(1-tz)*v1.x;
    v1.y = tz*v5.y+(1-tz)*v1.y;
    v1.z = tz*v5.z+(1-tz)*v1.z;

    P.x=v1.x;  P.y=v1.y;  P.z=v1.z;

    return 1;
}


//-------------------------------------------------------------

//----
void TensorField::createInterleavedDTI(ScalarField *Dxx, ScalarField *Dyy, ScalarField *Dzz, ScalarField *Dxy, ScalarField *Dxz, ScalarField *Dyz)
//----
{

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix, iy, iz;
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *Dxxval= Dxx->getAt(ix,iy,iz);
	  float *Dyyval= Dyy->getAt(ix,iy,iz);
	  float *Dzzval= Dzz->getAt(ix,iy,iz);
	  float *Dxyval= Dxy->getAt(ix,iy,iz);
	  float *Dxzval= Dxz->getAt(ix,iy,iz);
	  float *Dyzval= Dyz->getAt(ix,iy,iz);
	  
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  m_Data[ind*6]= Dxxval[0];
	  m_Data[ind*6+1]= Dyyval[0];
	  m_Data[ind*6+2]= Dzzval[0];
	  m_Data[ind*6+3]= Dxyval[0];
	  m_Data[ind*6+4]= Dxzval[0];
	  m_Data[ind*6+5]= Dyzval[0];

	}
}

//----
void TensorField::extractDTIcomponents(ScalarField *Dxx, ScalarField *Dyy, ScalarField *Dzz, ScalarField *Dxy, ScalarField *Dxz, ScalarField *Dyz)
//----
{

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  int ix, iy, iz;
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *Dxxval= Dxx->getAt(ix,iy,iz);
	  float *Dyyval= Dyy->getAt(ix,iy,iz);
	  float *Dzzval= Dzz->getAt(ix,iy,iz);
	  float *Dxyval= Dxy->getAt(ix,iy,iz);
	  float *Dxzval= Dxz->getAt(ix,iy,iz);
	  float *Dyzval= Dyz->getAt(ix,iy,iz);
	  
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  Dxxval[0]= m_Data[ind*6];
	  Dyyval[0]= m_Data[ind*6+1];
	  Dzzval[0]= m_Data[ind*6+2] ;
	  Dxyval[0]= m_Data[ind*6+3];
	  Dxzval[0]= m_Data[ind*6+4];
	  Dyzval[0]= m_Data[ind*6+5];

	}
}



//----
void TensorField::ComputeFA(ScalarField *faField)
//----
{
int ix, iy, iz, i, j;

  double M[9];
  double eigVal[3], fa, DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  gsl_vector *eval ;
  gsl_matrix *evec ;
  
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
   
  w = gsl_eigen_symmv_alloc (3);
    
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
     
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
    
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	  
	  DTtrace=0;
	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      DTtrace += eigVal[i];
	    }

	  
#ifdef _V2_
	  double DTtrace2= M[0]+M[4]+M[8];  
	  if(fabs(DTtrace2)<1e-20)
	    DTtrace= 0;
#endif

	  if(DTtrace>0)
	    fa = sqrt(3.0/2.0*((eigVal[0]-DTtrace/3)*(eigVal[0]-DTtrace/3) +  (eigVal[1]-DTtrace/3)*(eigVal[1]-DTtrace/3) + (eigVal[2]-DTtrace/3)*(eigVal[2]-DTtrace/3))/((eigVal[0]*eigVal[0]) +  (eigVal[1]*eigVal[1]) + (eigVal[2]*eigVal[2])));	  	  
	  else
	    fa=0;

	  // this can only happen due to negative eigenvalues
	  if(fa>1.0)
	    fa=1.0;

	  //faField->setVal(ind,fa);
	  float *val= faField->getAt(ix,iy,iz);
	  *val = fa;

	}

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  
}

//----
void TensorField::ComputeTrace(ScalarField *traceField)
//----
{
int ix, iy, iz, i, j;

  double  DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	
	  DTtrace=m_Data[ind*6] + m_Data[ind*6+1]+ m_Data[ind*6+2];
	  	 
	  float *val= traceField->getAt(ix,iy,iz);
	  *val = DTtrace;
	}
 
}

//----
void TensorField::ComputePD(VectorField *PDField, int d_index, int weight_flag)
//----
{
int ix, iy, iz, i, j;

  double M[9];
  double eigVal[3], fa, PD[3], DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];


  gsl_vector *eval ;
  gsl_matrix *evec ;
  
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
   
  w = gsl_eigen_symmv_alloc (3);
    
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
     
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
    
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	  
	  DTtrace=0;
	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);	      
	      DTtrace += eigVal[i];
	    }
	 
#ifdef _V2_
	  double DTtrace2= M[0]+M[4]+M[8];  
	  if(fabs(DTtrace2)<1e-20)
	    DTtrace= 0;
#endif


	  if(DTtrace>0)
	    {
	      /*
	       for(i=0; i<3; i++)
		 printf("%1.2f ",eigVal[i]);
	       printf("\n");
	      */

	      gsl_vector_view evec_i = gsl_matrix_column (evec, d_index);
	      
	      float PDnorm=0;
	      for(j=0; j<3; j++)
		{
		  PD[j] = gsl_vector_get (&evec_i.vector, j);		     
		  PDnorm += PD[j]*PD[j];
		}
	      
	      //std::cout<<"PDnorm: "<<PDnorm<<"\n"; //checked and always 1
	      if(weight_flag>0)	
		{
		  PDnorm=0.;
		   
		  fa = sqrt(3.0/2.0*((eigVal[0]-DTtrace/3)*(eigVal[0]-DTtrace/3) +  (eigVal[1]-DTtrace/3)*(eigVal[1]-DTtrace/3) + (eigVal[2]-DTtrace/3)*(eigVal[2]-DTtrace/3))/((eigVal[0]*eigVal[0]) +  (eigVal[1]*eigVal[1]) + (eigVal[2]*eigVal[2])));
	  	  
		  //clip fa: this can only be caused by negative eigenvalues
		  if(fa>1.0)
		    fa=1.0;
		  
		  for(j=0; j<3; j++)
		    {
		      PD[j] = fa*PD[j];
		      //   PDnorm += PD[j]*PD[j];
		    }
		  //  PDnorm = sqrt(PDnorm);
		  
		  //   std::cout<<"FA: "<<fa<<" "<<"PDnorm: "<<PDnorm<<"\n"; //checked and always equal
		}
	     
	    }
      
	  float *val= PDField->getAt(ix,iy,iz);
	  for(j=0; j<3; j++)
	    {
	      //PDField->setVal(ind*3+j,PD[j]);
	     
	      if(DTtrace>0)
		val[j] = PD[j];	     
	      else 
		val[j]= 0;
	    }
	  
	}
  
  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  
}
  
//----
void TensorField::ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field)
//----
{
  int ix, iy, iz, i, j;

  double M[9];
  double eigVal[3], eigVec[9], fa, DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];


  gsl_vector *eval ;
  gsl_matrix *evec ;
  
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
   
  w = gsl_eigen_symmv_alloc (3);
    
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
     
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
    
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	 
	  DTtrace= 0;
	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      DTtrace += eigVal[i];

	      gsl_vector_view evec_i 
                = gsl_matrix_column (evec, i);

	      for(j=0; j<3; j++)
		eigVec[i*3+j] = gsl_vector_get (&evec_i.vector, j);	     
		  
	    }

#ifdef _V2_
	  double DTtrace2= M[0]+M[4]+M[8];  
	  if(fabs(DTtrace2)<1e-20)
	    DTtrace= 0;
#endif


	  float *val;
	  //e1Field->setVal(ind, eigVal[0]);
	  val= e1Field->getAt(ix,iy,iz);
	  *val = eigVal[0];
	  //e2Field->setVal(ind, eigVal[1]);
	  val= e2Field->getAt(ix,iy,iz);
	  *val = eigVal[1];
	  //e3Field->setVal(ind, eigVal[2]);
	  val= e3Field->getAt(ix,iy,iz);
	  *val = eigVal[2];

	  for(j=0; j<3; j++)
	    {
	      //PD1Field->setVal(ind*3+j,eigVec[0*3+j]);
	      val= PD1Field->getAt(ix,iy,iz);
	      val[j] = eigVec[0*3+j];
	      //PD2Field->setVal(ind*3+j,eigVec[1*3+j]);
	      val= PD2Field->getAt(ix,iy,iz);
	      val[j] = eigVec[1*3+j];
	      //PD3Field->setVal(ind*3+j,eigVec[2*3+j]);
	      val= PD3Field->getAt(ix,iy,iz);
	      val[j] = eigVec[2*3+j];
	    }

#ifdef _V2_
	  if(DTtrace==0)
	    {
	       //e1Field->setVal(ind, eigVal[0]);
	      val= e1Field->getAt(ix,iy,iz);
	      *val = 0;
	      //e2Field->setVal(ind, eigVal[1]);
	      val= e2Field->getAt(ix,iy,iz);
	      *val = 0;
	      //e3Field->setVal(ind, eigVal[2]);
	      val= e3Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      for(j=0; j<3; j++)
		{
		  //PD1Field->setVal(ind*3+j,eigVec[0*3+j]);
		  val= PD1Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD2Field->setVal(ind*3+j,eigVec[1*3+j]);
		  val= PD2Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD3Field->setVal(ind*3+j,eigVec[2*3+j]);
		  val= PD3Field->getAt(ix,iy,iz);
		  val[j] = 0;
		}
	    }  
#endif

	}	 

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  
}



//----
void TensorField::ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field, ScalarField *traceField, ScalarField *faField, int FAweighting)
//----
{
  int ix, iy, iz, i, j;

  double M[9];
  double eigVal[3], eigVec[9], fa, DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];


  gsl_vector *eval ;
  gsl_matrix *evec ;
  
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
   
  w = gsl_eigen_symmv_alloc (3);
    
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
     
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
    
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	 
	  DTtrace= 0;
	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      DTtrace += eigVal[i];

	      gsl_vector_view evec_i 
                = gsl_matrix_column (evec, i);

	      for(j=0; j<3; j++)
		eigVec[i*3+j] = gsl_vector_get (&evec_i.vector, j);	     
		  
	    }

#ifdef _V2_
	  double DTtrace2= M[0]+M[4]+M[8];  
	  if(fabs(DTtrace2)<1e-20)
	    DTtrace= 0;
#endif


	  float *val;
	  
	  val= e1Field->getAt(ix,iy,iz);
	  *val = eigVal[0];
	  
	  val= e2Field->getAt(ix,iy,iz);
	  *val = eigVal[1];
	  
	  val= e3Field->getAt(ix,iy,iz);
	  *val = eigVal[2];

	  val= traceField->getAt(ix,iy,iz);
	  *val = DTtrace;

	  if(DTtrace>0)
	    fa = sqrt(3.0/2.0*((eigVal[0]-DTtrace/3)*(eigVal[0]-DTtrace/3) +  (eigVal[1]-DTtrace/3)*(eigVal[1]-DTtrace/3) + (eigVal[2]-DTtrace/3)*(eigVal[2]-DTtrace/3))/((eigVal[0]*eigVal[0]) +  (eigVal[1]*eigVal[1]) + (eigVal[2]*eigVal[2])));	  	  
	  else
	    fa=0;

	  // this can only happen due to negative eigenvalues
	  if(fa>1.0)
	    fa=1.0;

	  //faField->setVal(ind,fa);
	  val= faField->getAt(ix,iy,iz);
	  *val = fa;

	  for(j=0; j<3; j++)
	    {
	      val= PD1Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[0*3+j];
	      else
		val[j] = eigVec[0*3+j];
	  
	      val= PD2Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[1*3+j];
	      else
		val[j] = eigVec[1*3+j];
	  
	      val= PD3Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[2*3+j];
	      else
		val[j] = eigVec[2*3+j];
	    }
	  

#ifdef _V2_
	  if(DTtrace==0)
	    {
	      
	      val= e1Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= e2Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= e3Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= traceField->getAt(ix,iy,iz);
	      *val = 0;

	      val= faField->getAt(ix,iy,iz);
	      *val = 0;

	      for(j=0; j<3; j++)
		{
		  //PD1Field->setVal(ind*3+j,eigVec[0*3+j]);
		  val= PD1Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD2Field->setVal(ind*3+j,eigVec[1*3+j]);
		  val= PD2Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD3Field->setVal(ind*3+j,eigVec[2*3+j]);
		  val= PD3Field->getAt(ix,iy,iz);
		  val[j] = 0;
		}
	    }  
#endif


	} 

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  
}


void TensorField::EulerAngles2RotMat(double *RotMat, float phi, float theta, float psi)
{

  RotMat[0*3+0]= cos(theta)*cos(phi);
  RotMat[1*3+0]= cos(theta)*sin(phi);
  RotMat[2*3+0]= -sin(theta);
  
  RotMat[0*3+1]= sin(psi)*sin(theta)*cos(phi)- cos(psi)*sin(phi);
  RotMat[1+3+1]= sin(psi)*sin(theta)*sin(phi)+ cos(psi)*cos(phi);
  RotMat[2*3+1]= sin(psi)*cos(theta);
  
  RotMat[0*3+2]= cos(psi)*sin(theta)*cos(phi)+ sin(psi)*sin(phi);
  RotMat[1*3+2]= cos(psi)*sin(theta)*sin(phi)- sin(psi)*cos(phi);
  RotMat[2*3+2]= cos(psi)*cos(theta);
 
}

void TensorField::RotMat2EulerAngles(const double *RotMat, float &phi, float &theta, float &psi)
{
  
  float phi1, theta1, psi1, phi2, theta2, psi2, delta;
  
  if(fabs(fabs(RotMat[2*3+0])-1)>1e-20)
    {
      theta1= -asin(RotMat[2*3+0]);
      theta2= M_PI - theta1;
      psi1= atan2(RotMat[2*3+1]/cos(theta1),RotMat[2*3+2]/cos(theta1));
      psi2= atan2(RotMat[2*3+1]/cos(theta2),RotMat[2*3+2]/cos(theta2));
      phi1= atan2(RotMat[1*3+0]/cos(theta1),RotMat[0*3+0]/cos(theta1));
      phi2= atan2(RotMat[1*3+0]/cos(theta2),RotMat[0*3+0]/cos(theta2));
      
      phi= phi1;
      theta= theta1;
      psi= psi1;      
    }
  else
    {
      phi = 0;
      delta = atan2(RotMat[0*3+1],RotMat[0*3+2]);
      if(fabs(RotMat[2*3+0]+1) < 1e-20)
	{
	  theta= M_PI/2;
	  psi= phi + delta;
	}
      else
	{
	  theta= -M_PI/2;
	  psi= -phi + delta;
	}
    }
}

//----
void TensorField::ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field, ScalarField *traceField, ScalarField *faField, int FAweighting, ScalarField *PhiField, ScalarField *ThetaField, ScalarField *PsiField)
//----
{
  int ix, iy, iz, i, j;

  double M[9];
  double eigVal[3], eigVec[9], fa, DTtrace;

  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];

  gsl_vector *eval ;
  gsl_matrix *evec ;
  
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
   
  w = gsl_eigen_symmv_alloc (3);
    
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
     
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
    
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	 
	  DTtrace= 0;
	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      DTtrace += eigVal[i];

	      gsl_vector_view evec_i 
                = gsl_matrix_column (evec, i);

	      for(j=0; j<3; j++)
		eigVec[i*3+j] = gsl_vector_get (&evec_i.vector, j);	     
		  
	    }

#ifdef _V2_
	  double DTtrace2= M[0]+M[4]+M[8];  
	  if(fabs(DTtrace2)<1e-20)
	    DTtrace= 0;
#endif



	  float *val;
	  
	  val= e1Field->getAt(ix,iy,iz);
	  *val = eigVal[0];
	  
	  val= e2Field->getAt(ix,iy,iz);
	  *val = eigVal[1];
	  
	  val= e3Field->getAt(ix,iy,iz);
	  *val = eigVal[2];

	  val= traceField->getAt(ix,iy,iz);
	  *val = DTtrace;

	  if(DTtrace>0)
	    fa = sqrt(3.0/2.0*((eigVal[0]-DTtrace/3)*(eigVal[0]-DTtrace/3) +  (eigVal[1]-DTtrace/3)*(eigVal[1]-DTtrace/3) + (eigVal[2]-DTtrace/3)*(eigVal[2]-DTtrace/3))/((eigVal[0]*eigVal[0]) +  (eigVal[1]*eigVal[1]) + (eigVal[2]*eigVal[2])));	  	  
	  else
	    fa=0;

	  // this can only happen due to negative eigenvalues
	  if(fa>1.0)
	    fa=1.0;

	  //faField->setVal(ind,fa);
	  val= faField->getAt(ix,iy,iz);
	  *val = fa;

	  for(j=0; j<3; j++)
	    {
	      val= PD1Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[0*3+j];
	      else
		val[j] = eigVec[0*3+j];
	  
	      val= PD2Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[1*3+j];
	      else
		val[j] = eigVec[1*3+j];
	  
	      val= PD3Field->getAt(ix,iy,iz);
	      if(FAweighting)
		val[j] = fa*eigVec[2*3+j];
	      else
		val[j] = eigVec[2*3+j];
	    }
	  
	  float phi, theta, psi;
	  RotMat2EulerAngles(eigVec, phi, theta, psi);
	  val= PhiField->getAt(ix,iy,iz);
	  val[0]= phi;
	  val= ThetaField->getAt(ix,iy,iz);
	  val[0]= theta;
	  val= PsiField->getAt(ix,iy,iz);
	  val[0]= psi;

#ifdef _V2_
	  if(DTtrace==0)
	    {
	      
	      val= e1Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= e2Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= e3Field->getAt(ix,iy,iz);
	      *val = 0;
	      
	      val= traceField->getAt(ix,iy,iz);
	      *val = 0;

	      val= faField->getAt(ix,iy,iz);
	      *val = 0;

	      for(j=0; j<3; j++)
		{
		  //PD1Field->setVal(ind*3+j,eigVec[0*3+j]);
		  val= PD1Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD2Field->setVal(ind*3+j,eigVec[1*3+j]);
		  val= PD2Field->getAt(ix,iy,iz);
		  val[j] = 0;
		  //PD3Field->setVal(ind*3+j,eigVec[2*3+j]);
		  val= PD3Field->getAt(ix,iy,iz);
		  val[j] = 0;
		}

	       val= PhiField->getAt(ix,iy,iz);
	       val[0]= 0;
	       val= ThetaField->getAt(ix,iy,iz);
	       val[0]= 0;
	       val= PsiField->getAt(ix,iy,iz);
	       val[0]= 0;


	    }  
#endif

	  

	} 

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  
}



//----
void TensorField::computeNZmask(unsigned char *mask, float threshold_factor)
//----
{
  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float *dti_data = getAt(Xdim/2,Ydim/2,Zdim/2);
  float cen_trace = dti_data[0]+dti_data[1]+dti_data[2];

  traceTensorThreshold = threshold_factor*cen_trace ;

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *val = getAt(ix,iy,iz);
	  float trace = val[0]+val[1]+val[2];
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  if(trace > traceTensorThreshold)
	    mask[ind]=1;
	  else
	    mask[ind]=0;
	}
}


//----
void TensorField::logTField(TensorField *dti_out, unsigned char *mask)
//----
{
  int ix, iy, iz, i, j;
  
  double M[9];
  double eigVal[3], eigVec[9], logMat[9], fa, DTtrace;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  /*
  TensorField dti_out;
  
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  */  
  
  gsl_vector *eval ;
  gsl_matrix *evec ;
  gsl_matrix *evec_trans ;
  gsl_matrix *eval_mat ;
  gsl_matrix *temp_mat ;
  gsl_matrix *log_mat ;
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
  evec_trans = gsl_matrix_alloc (3, 3);
  eval_mat = gsl_matrix_alloc(3,3);
  temp_mat = gsl_matrix_alloc(3,3);
  log_mat = gsl_matrix_alloc(3,3);  

  w = gsl_eigen_symmv_alloc (3);
  
  // first compute log of tensor 
  // we store it temporarily into log
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];	  	  

	  DTtrace= M[0] +  M[4] +  M[8];

	  if(mask[ind]==0 || DTtrace < traceTensorThreshold)
	    {
	      float *val= dti_out->getAt(ix,iy,iz);
	      for(j=0; j<6; j++)
		val[j]=0;
	      continue;
	    }


	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
	  
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
	  
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	  
	 
	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      gsl_matrix_set(eval_mat, i,j,0.0);


	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      
	      if(eigVal[i]>0.0)		
		gsl_matrix_set(eval_mat, i,i,log(eigVal[i]));

	      /*--
	      gsl_vector_view evec_i 
		= gsl_matrix_column (evec, i);
	      
	      for(j=0; j<3; j++)
		eigVec[i*3+j] = gsl_vector_get (&evec_i.vector, j);	     
		--*/
	    }
	  
	  //form the log tensor and insert into dti_out
	  //evec*eval_mat*evec_t;
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, eval_mat, 0.0, temp_mat);	  
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_mat, evec, 0.0, log_mat);
	  for(i=0; i<3; i++)
	    {	     	     
	      for(j=0; j<3; j++)
		logMat[i*3+j] = gsl_matrix_get (log_mat, i, j);	     
	    }
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  val[0]= (float) logMat[0] ;
	  val[1]= (float) logMat[4] ;
	  val[2]= (float) logMat[8] ;
	  val[3]= (float) logMat[1] ;
	  val[4]= (float) logMat[2] ;
	  val[5]= (float) logMat[5] ;

	}
 

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  gsl_matrix_free( eval_mat);
  gsl_matrix_free( temp_mat);
  gsl_matrix_free( log_mat);

  //  return(dti_out);

}

//----
void TensorField::expTField(TensorField *dti_out, unsigned char *mask)
//----
{
  int ix, iy, iz, i, j;
  
  double M[9];
  double eigVal[3], eigVec[9], expMat[9], fa;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  /*
  TensorField dti_out;
  
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  */
  
  gsl_vector *eval ;
  gsl_matrix *evec ;
  gsl_matrix *evec_trans ;
  gsl_matrix *eval_mat ;
  gsl_matrix *temp_mat ;
  gsl_matrix *exp_mat ;
  gsl_eigen_symmv_workspace  *w ;
  
  eval = gsl_vector_alloc (3);
  evec = gsl_matrix_alloc (3, 3);
  evec_trans = gsl_matrix_alloc (3, 3);
  eval_mat = gsl_matrix_alloc(3,3);
  temp_mat = gsl_matrix_alloc(3,3);
  exp_mat = gsl_matrix_alloc(3,3);  

  w = gsl_eigen_symmv_alloc (3);
  
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];
	  
	  if(mask[ind]==0 || (M[0] +  M[4] +  M[8])==0)
	    {
	      float *val= dti_out->getAt(ix,iy,iz);
	      for(j=0; j<6; j++)
		val[j]=0;
	      continue;
	    }


	  gsl_matrix_view m 
	    = gsl_matrix_view_array (M, 3, 3);
	  
	  gsl_eigen_symmv (&m.matrix, eval, evec, w);
	  
	  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	  
	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      gsl_matrix_set(eval_mat, i,j,0.0);

	  for(i=0; i<3; i++)
	    {
	      eigVal[i] = gsl_vector_get (eval, i);
	      
	      gsl_matrix_set(eval_mat, i,i,exp(eigVal[i]));

	      /*--
	      gsl_vector_view evec_i 
		= gsl_matrix_column (evec, i);
	      
	      for(j=0; j<3; j++)
		eigVec[i*3+j] = gsl_vector_get (&evec_i.vector, j);	     
		--*/
	    }
	  
	  //form the exp tensor and insert into dti_out
	  //evec*eval_mat*evec_t;
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, eval_mat, 0.0, temp_mat);	  
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_mat, evec, 0.0, exp_mat);
	  for(i=0; i<3; i++)
	    {	     	     
	      for(j=0; j<3; j++)
		expMat[i*3+j] = gsl_matrix_get (exp_mat, i, j);	     
	    }
	  	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  val[0]= (float) expMat[0] ;
	  val[1]= (float) expMat[4] ;
	  val[2]= (float) expMat[8] ;
	  val[3]= (float) expMat[1] ;
	  val[4]= (float) expMat[2] ;
	  val[5]= (float) expMat[5] ;

	}

  gsl_eigen_symmv_free ( w);
  gsl_vector_free( eval);
  gsl_matrix_free( evec);
  gsl_matrix_free( eval_mat);
  gsl_matrix_free( temp_mat);
  gsl_matrix_free( exp_mat);

  // return(dti_out);

}

//----
void TensorField::ESmooth(ScalarField *smoothingFilter, TensorField *dti_out)
//----
{

  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  

  int ix1, iy1, iz1;
  int ix2, iy2, iz2;

  int fXdim = smoothingFilter->getSize(0);
  int fYdim = smoothingFilter->getSize(1);
  int fZdim = smoothingFilter->getSize(2);

  /*--
  TensorField dti_out;
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  --*/

  if(fXdim%2 ==0) { std::cout<<"Kernel size should be odd\n"; exit(0); }
  if(fYdim%2 ==0) { std::cout<<"Kernel size should be odd\n"; exit(0); }
  if(fZdim%2 ==0) { std::cout<<"Kernel size should be odd\n"; exit(0); }

  // assume kernel sizes are odd
  
  float *kernel_Data = smoothingFilter->getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  for(ix=0; ix<Xdim; ix++)
    {
      for(iy=0; iy<Ydim; iy++)
	for(iz=0; iz<Zdim; iz++)
	  {
	  
	 
	  float *val= dti_out->getAt(ix,iy,iz);
	  for(j=0; j<6; j++)
	    val[j] = 0.0;

	  if(ix < order_x || ix >= Xdim- order_x || iy<order_y || iy >= Ydim- order_y || iz<order_z || iz>=Zdim- order_z )
	    continue;

	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;

		  ix2 = ix -ix1 + order_x;
		  iy2 = iy -iy1 + order_y;
		  iz2 = iz -iz1 + order_z;
		  
		  int ind = iz2*Xdim*Ydim + iy2*Xdim + ix2;
		  
		  for(j=0; j<6; j++)
		    val[j] += kernel_Data[ind1]*m_Data[ind*6+j];
     		}
	  }
      //std::cout<<ix<<" YZ slice smoothing done\n";
    }
  
}

//----
void TensorField::ESmooth(ScalarField *smoothingFilter, TensorField *dti_out,  unsigned char *mask)
//----
{

  int ix, iy, iz, i, j;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  int ix1, iy1, iz1;
  int ix2, iy2, iz2;

  int fXdim = smoothingFilter->getSize(0);
  int fYdim = smoothingFilter->getSize(1);
  int fZdim = smoothingFilter->getSize(2);

  

  /*--
  TensorField dti_out;
  dti_out.init(Xdim,Ydim,Zdim);
  dti_out.setVoxelSize(0,Xres); dti_out.setVoxelSize(1,Yres); dti_out.setVoxelSize(2,Zres);
  --*/

  if(fXdim%2 ==0) { std::cout<<"ESmooth: Kernel size should be odd\n"; exit(0); }
  if(fYdim%2 ==0) { std::cout<<"ESmooth: Kernel size should be odd\n"; exit(0); }
  if(fZdim%2 ==0) { std::cout<<"ESmooth: Kernel size should be odd\n"; exit(0); }


  // first smooth the mask
  ScalarField mask_float;
  mask_float.init(Xdim,Ydim,Zdim);
  mask_float.setVoxelSize(0,Xres); 
  mask_float.setVoxelSize(1,Yres); 
  mask_float.setVoxelSize(2,Zres);

  ScalarField mask_float_smooth;
  mask_float_smooth.init(Xdim,Ydim,Zdim);
  mask_float_smooth.setVoxelSize(0,Xres); 
  mask_float_smooth.setVoxelSize(1,Yres); 
  mask_float_smooth.setVoxelSize(2,Zres);


  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *mask_val = mask_float.getAt(ix,iy,iz);
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  mask_val[0] = (float) mask[ind];
	}
	
  mask_float.Smooth(smoothingFilter, &mask_float_smooth);
  // std::cout<<"Smooth mask computed\n";

  // assume kernel sizes are odd
  
  float *kernel_Data = smoothingFilter->getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  for(ix=0; ix<Xdim; ix++)
    {
      for(iy=0; iy<Ydim; iy++)
	for(iz=0; iz<Zdim; iz++)
	  {
	  
	 
	  float *val= dti_out->getAt(ix,iy,iz);
	  for(j=0; j<6; j++)
	    val[j] = 0.0;

	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  if(mask[ind]==0)
	    continue;

	  if(ix < order_x || ix >= Xdim- order_x || iy<order_y || iy >= Ydim- order_y || iz<order_z || iz>=Zdim- order_z )
	    continue;

	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;

		  ix2 = ix -ix1 + order_x;
		  iy2 = iy -iy1 + order_y;
		  iz2 = iz -iz1 + order_z;
		  
		  int ind = iz2*Xdim*Ydim + iy2*Xdim + ix2;
		  
		  for(j=0; j<6; j++)
		    val[j] += kernel_Data[ind1]*m_Data[ind*6+j];
     		}
	  }
      //std::cout<<ix<<" YZ slice smoothing done\n";
    }
  
  // normalize with smooth mask
  for(ix=0; ix<Xdim; ix++)
    {
      for(iy=0; iy<Ydim; iy++)
	for(iz=0; iz<Zdim; iz++)
	  {
	    float *val= dti_out->getAt(ix,iy,iz);
	    
	    int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	    if(mask[ind]==0)
	      continue;

	    float *mask_smooth_val = mask_float_smooth.getAt(ix,iy,iz);
	    
	    if(mask_smooth_val[0]>0)
	      for(j=0; j<6; j++)
		val[j] /= mask_smooth_val[0];
	  }
    }
  // std::cout<<"Normalization with smooth mask done\n";

}

//----
void TensorField::logESmooth(ScalarField *smoothingFilter, TensorField *dti_out)
//----
{
  int ix, iy, iz, i, j;
  
  double M[9];
  double eigVal[3], eigVec[9], logMat[9], fa;
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  TensorField dti_out1, dti_out2;
  
  dti_out1.init(Xdim,Ydim,Zdim);
  dti_out1.setVoxelSize(0,Xres); 
  dti_out1.setVoxelSize(1,Yres); 
  dti_out1.setVoxelSize(2,Zres);
  
  dti_out2.init(Xdim,Ydim,Zdim);
  dti_out2.setVoxelSize(0,Xres); 
  dti_out2.setVoxelSize(1,Yres); 
  dti_out2.setVoxelSize(2,Zres);

  unsigned char *mask= (unsigned char *)calloc(Xdim*Ydim*Zdim, sizeof(unsigned char));

  // note magic number 0.001
  computeNZmask(mask,0.001);

  char str[200];

  // first compute log of tensor 
  logTField(&dti_out1, mask);
  //std::cout<<"Log Tensor Field computed\n";  

  /*--
  sprintf(str,"/sbia/home/pkhurd/XDR2HariDTcode/DATA/logDTI%d.img",Xdim);
  dti_out1.exportRawFieldToFile(str);

  dti_out1.expTField(dti_out, mask ); // should recover original
  sprintf(str,"/sbia/home/pkhurd/XDR2HariDTcode/DATA/origDTI%d.img",Xdim);
  dti_out->exportRawFieldToFile(str);
  --*/

  // now do Euclidean smoothing
  dti_out1.ESmooth(smoothingFilter, &dti_out2, mask);
  //std::cout<<"Log Tensor Field smoothing done\n";  

  // now exponentiate tensor
  //dti_out2.computeNZmask(mask,0.1);
  dti_out2.expTField(dti_out, mask );
  //std::cout<<"Exp Tensor Field computed\n";  

  free(mask);

}

//----
void TensorField::reOrientTensorFieldFS(VectorField *warpField, TensorField *dti_out)
//----
{
  int ix, iy, iz, i, j;
 
  double M[9];
  
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
     
  double U[9]; double outMat[9];
  
  gsl_vector *temp_vec1 ;
  gsl_vector *temp_vec2 ;
  gsl_matrix *temp_mat1 ;
  gsl_matrix *temp_mat2 ;
  gsl_matrix *temp_mat3 ;
  gsl_eigen_symmv_workspace  *w ;
  
  temp_vec1 = gsl_vector_alloc(3);
  temp_vec2 = gsl_vector_alloc(3);
  temp_mat1 = gsl_matrix_alloc(3,3);
  temp_mat2 = gsl_matrix_alloc(3,3);
  temp_mat3 = gsl_matrix_alloc(3,3);  
  w = gsl_eigen_symmv_alloc (3);
  
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	
	  //Step 1: estimate rotation with the Finite Strain approximation----
	  
	  // Step 1a: get Jacobian
	  // here we assume warp field vectors map old grid to new grid
	  // here we store transpose of Jacobian in temp_mat1
	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      {
		if(i==j)
		  gsl_matrix_set(temp_mat1,i,j,1.0);
		else
		  gsl_matrix_set(temp_mat1,i,j,0.0);
	      }

	  if(ix<Xdim-1 && iy<Ydim-1 && iz<Zdim-1)
	    {
	      float *warpField_val = warpField->getAt(ix,iy,iz);
	      float *warpField_val2;
	      for(i=0; i<3; i++) // index column: old grid
		{ 
		  switch(i)
		    {
		    case 0:
		      warpField_val2 = warpField->getAt(ix+1,iy,iz);
		      break;
		    case 1:
		      warpField_val2 = warpField->getAt(ix,iy+1,iz);
		      break;
		    case 2:
		      warpField_val2 = warpField->getAt(ix,iy,iz+1);
		      break;
		    }
		  for(j=0; j<3; j++) // index row: new grid
		    {
		      float diff = warpField_val2[j] - warpField_val[j]; 
		      gsl_matrix_set(temp_mat1,i,j,diff);
		    }

		   if(ix== Xdim/2+13 && iy== Ydim/2+13 && iz==Zdim/2+13)
		     {
		       std::cout<<"Field at "<<ix << "," << iy << "," <<iz<< "\n";
		       for(j=0; j<3; j++)
			 std::cout<<warpField_val2[j]<<" ";
		       std::cout<<"\n";
		     }
		}
	    }
	    
	  // for debug
	  if(ix== Xdim/2+13 && iy== Ydim/2+13 && iz==Zdim/2+13)
	    {
	      
	      std::cout<<"Estimated Jacobian at "<< ix << "," << iy << "," <<iz<< "\n";
	      for(int pi=0; pi<3; pi++)
		std::cout<<gsl_matrix_get(temp_mat1,0,pi)<<" "<<gsl_matrix_get(temp_mat1,1,pi)<<" "<<gsl_matrix_get(temp_mat1,2,pi)<<"\n";
	      
	    }

	  // Step 1b: get rotation component

	  gsl_linalg_SV_decomp (temp_mat1, temp_mat3, temp_vec1, temp_vec2);
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_mat1, temp_mat3, 0.0, temp_mat2);

	    /*-- this was a wrong idea --- QR of Jac*OrigEV is related to PPD 
	  gsl_linalg_QR_decomp (temp_mat1, temp_vec1);
	  gsl_linalg_QR_unpack (temp_mat1, temp_vec1, temp_mat2, temp_mat3);
	  --*/

	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      U[i*3+j] = gsl_matrix_get(temp_mat2, j,i);
	      //U[i*3+j] = gsl_matrix_get(temp_mat2, i,j);
	  

	  // for debug
	  if(ix== Xdim/2+13 && iy== Ydim/2+13 && iz==Zdim/2+13)
	    {
	      std::cout<<"Estimated Rotation at "<<ix<<","<<iy<<","<<iz<<"\n";
	      for(int pi=0; pi<3; pi++)
		std::cout<<gsl_matrix_get(temp_mat2,0,pi)<<" "<<gsl_matrix_get(temp_mat2,1,pi)<<" "<<gsl_matrix_get(temp_mat2,2,pi)<<"\n";
	      
	    }


	  //Step 2: apply rotation ----

	  M[0]= (double) m_Data[ind*6];
	  M[4]= (double) m_Data[ind*6+1];
	  M[8]= (double) m_Data[ind*6+2];
	  M[1]= (double) m_Data[ind*6+3];
	  M[2]= (double) m_Data[ind*6+4];
	  M[5]= (double) m_Data[ind*6+5];
	  M[3]= M[1];
	  M[6]= M[2];
	  M[7]= M[5];

	   for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      {
		gsl_matrix_set(temp_mat1, i,j, U[i*3+j]); //rotation
		gsl_matrix_set(temp_mat2, i,j, M[i*3+j]); // original tensor
	      }

	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp_mat1, temp_mat2, 0.0, temp_mat3);
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_mat3, temp_mat1, 0.0, temp_mat2);

	  
	  for(i=0; i<3; i++)
	    {	     	     
	      for(j=0; j<3; j++)
		outMat[i*3+j] = gsl_matrix_get (temp_mat2, i, j);	     
	    }
	  
	   // for debug
	  if(ix== Xdim/2+13 && iy== Ydim/2+13 && iz==Zdim/2+13)
	    {
	      std::cout<<"Original Tensor at "<<ix<<","<<iy<<","<<iz<<"\n";
	      for(int pi=0; pi<3; pi++)
		std::cout<<M[pi*3]<<" "<<M[pi*3+1]<<" "<<M[pi*3+2]<<"\n";

	      std::cout<<"Rotated Tensor at "<<ix<<","<<iy<<","<<iz<<"\n";
	      for(int pi=0; pi<3; pi++)
		std::cout<<outMat[pi*3]<<" "<<outMat[pi*3+1]<<" "<<outMat[pi*3+2]<<"\n";
	      
	    }

	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  val[0]= (float) outMat[0] ;
	  val[1]= (float) outMat[4] ;
	  val[2]= (float) outMat[8] ;
	  val[3]= (float) outMat[1] ;
	  val[4]= (float) outMat[2] ;
	  val[5]= (float) outMat[5] ;


	}

  gsl_eigen_symmv_free ( w);  
  gsl_vector_free( temp_vec1);
  gsl_vector_free( temp_vec2);
  gsl_matrix_free( temp_mat1);
  gsl_matrix_free( temp_mat2);
  gsl_matrix_free( temp_mat3);

}

//----
void TensorField::displaceTensorFieldUsingInvFieldTL(VectorField *warpField, TensorField *dti_out)
//----
// here, we warp subject to template using template2subject field
//  we use tri-linear interpolation here
{
  int ix, iy, iz, i, j;
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];

  
  int ix_subject, iy_subject, iz_subject;
  int ix1, iy1, iz1, ix2, iy2, iz2;
 
  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  for(j=0; j<6; j++)
	    val[j] = 0.0f;

	 
	  float *warpField_val = warpField->getAt(ix,iy,iz);
	  ix_subject = floor(warpField_val[0]);
	  iy_subject = floor(warpField_val[1]);
	  iz_subject = floor(warpField_val[2]);
	  
	  if(ix_subject < 0 || ix_subject >= Xdim || iy_subject< 0 || iy_subject >= Ydim || iz_subject< 0 || iz_subject>=Zdim)
	    {	      
	      continue;
	    }
	
	  float tl_x = warpField_val[0]- ix_subject;
	  float tl_y = warpField_val[1]- iy_subject;
	  float tl_z = warpField_val[2]- iz_subject;
				
	  float ww[8];
	  
	  ww[0] = (1-tl_x)*(1-tl_y)*(1-tl_z);
	  ww[1] = tl_x*(1-tl_y)*(1-tl_z);
	  ww[2] = tl_x*tl_y*(1-tl_z);
	  ww[3] = (1-tl_x)*tl_y*(1-tl_z);
	  ww[4] = (1-tl_x)*(1-tl_y)*tl_z;
	  ww[5] = tl_x*(1-tl_y)*tl_z;
	  ww[6] = tl_x*tl_y*tl_z;
	  ww[7] = (1-tl_x)*tl_y*tl_z;
	  
	  for(i=0; i<6; i++)		  
	    {
	      float *dt_val0= this->getAt(ix_subject, iy_subject, iz_subject);
	      float *dt_val1= this->getAt(ix_subject+1, iy_subject, iz_subject);
	      float *dt_val2= this->getAt(ix_subject+1, iy_subject+1, iz_subject);
	      float *dt_val3= this->getAt(ix_subject, iy_subject+1, iz_subject);
	      float *dt_val4= this->getAt(ix_subject, iy_subject, iz_subject+1);
	      float *dt_val5= this->getAt(ix_subject+1, iy_subject, iz_subject+1);
	      float *dt_val6= this->getAt(ix_subject+1, iy_subject+1, iz_subject+1);
	      float *dt_val7= this->getAt(ix_subject, iy_subject+1, iz_subject+1);
	      
	      val[i] = ww[0]*dt_val0[i] + ww[1]*dt_val1[i] + ww[2]*dt_val2[i] + ww[3]*dt_val3[i] + ww[4]*dt_val4[i] + ww[5]*dt_val5[i] + ww[6]*dt_val6[i] + ww[7]*dt_val7[i] ;

	    }
  	
	}

}


//----
void TensorField::displaceTensorFieldUsingInvField(VectorField *warpField, TensorField *dti_out, float sigma, float filter_size_red_factor)
//----
// here, we warp subject to template using template2subject field
//  do trunc-Gaussian interpolation with appropriate normalization
// (in future we may do trunc-Gaussian interpolation in Log-Euclidean domain) 
{
  int ix, iy, iz, i, j;
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];

  int fXdim = Xdim/filter_size_red_factor/Xres;
  int fYdim = Ydim/filter_size_red_factor/Yres;
  int fZdim = Zdim/filter_size_red_factor/Zres;
   
  if(fXdim%2==0) fXdim +=1;
  if(fYdim%2==0) fYdim +=1;
  if(fZdim%2==0) fZdim +=1;
  
  ScalarField smoothFilter;
  
  smoothFilter.init(fXdim,fYdim,fZdim);
  smoothFilter.setVoxelSize(0,Xres); 
  smoothFilter.setVoxelSize(1,Yres); 
  smoothFilter.setVoxelSize(2,Zres);
 
  smoothFilter.ComputeSmoothingFilter(sigma);
    
  // assume kernel sizes are odd
  
  float *kernel_Data = smoothFilter.getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  int ix_subject, iy_subject, iz_subject;
  int ix1, iy1, iz1, ix2, iy2, iz2;

  // create a unity mask for normalization
  
  ScalarField mask_float;
  mask_float.init(Xdim,Ydim,Zdim);
  mask_float.setVoxelSize(0,Xres); 
  mask_float.setVoxelSize(1,Yres); 
  mask_float.setVoxelSize(2,Zres);

  ScalarField mask_float_smooth;
  mask_float_smooth.init(Xdim,Ydim,Zdim);
  mask_float_smooth.setVoxelSize(0,Xres); 
  mask_float_smooth.setVoxelSize(1,Yres); 
  mask_float_smooth.setVoxelSize(2,Zres);


  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *mask_val = mask_float.getAt(ix,iy,iz);	  
	  mask_val[0] = (float) 1.0;
	}
	
  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  for(j=0; j<6; j++)
	    val[j] = 0.0f;

	  float *mask_float_smooth_val = mask_float_smooth.getAt(ix,iy,iz);
	  
	  float *warpField_val = warpField->getAt(ix,iy,iz);
	  ix_subject = rint(warpField_val[0]);
	  iy_subject = rint(warpField_val[1]);
	  iz_subject = rint(warpField_val[2]);

	  
	  if(ix_subject < 0 || ix_subject >= Xdim || iy_subject< 0 || iy_subject >= Ydim || iz_subject< 0 || iz_subject>=Zdim)
	    {
	      //std::cout<<"Voxel ("<<ix<<","<<iy<<","<<iz<<") mapped outside\n";
	      continue;
	    }
	  
	  smoothFilter.ComputeSmoothingFilter(sigma, warpField_val[0]-ix_subject + order_x, warpField_val[1]-iy_subject + order_y, warpField_val[2]-iz_subject + order_z );    


	  // do trunc-Gaussian interpolation
	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;
		  
		  ix2 = ix_subject -ix1 + order_x;
		  iy2 = iy_subject -iy1 + order_y;
		  iz2 = iz_subject -iz1 + order_z;
		  
		  if(ix2 < 0 || ix2 >= Xdim || iy2< 0 || iy2 >= Ydim || iz2< 0 || iz2>=Zdim)
		    continue;

		  float *val_subject= getAt(ix2,iy2,iz2);
		  
		  for(j=0; j<6; j++)
		    val[j] += kernel_Data[ind1]*val_subject[j];


		  float *mask_float_val = mask_float.getAt(ix2,iy2,iz2);

		  mask_float_smooth_val[0] += kernel_Data[ind1]*mask_float_val[0];
		}
	  
	  // pointwise division for normalization
	  if(mask_float_smooth_val[0] > 0)
	    for(j=0; j<6; j++)
	      val[j] /= mask_float_smooth_val[0];

	}

  // drop all tensors where trace falls below a threshold

  float *test_trace_val = dti_out->getAt(Zdim/2,Ydim/2,Zdim/2);
  float trace_cen = test_trace_val[0] + test_trace_val[1] + test_trace_val[2] ;

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  
	  float trace = val[0]+val[1]+val[2];

	  if(trace < 0.01*trace_cen)
	    
	  for(j=0; j<6; j++)
	    val[j] = 0.0f;

	}


}


//----
void TensorField::displaceTensorFieldUsingFwdField(VectorField *warpField, TensorField *dti_out, float sigma, float filter_size_red_factor)
//----
// here, we warp subject to template using subject2template field
// do trunc-Gaussian interpolation  with appropriate normalization
// (in future we may do trunc-Gaussian interpolation in Log-Euclidean domain)  
{
  int ix, iy, iz, i, j;
  int Xdim= m_Dims[0];
  int Ydim= m_Dims[1];
  int Zdim= m_Dims[2];
  
  float Xres = m_VoxelSize[0];
  float Yres = m_VoxelSize[1];
  float Zres = m_VoxelSize[2];
  
  int fXdim = Xdim/filter_size_red_factor/Xres;
  int fYdim = Ydim/filter_size_red_factor/Yres;
  int fZdim = Zdim/filter_size_red_factor/Zres;
   
  if(fXdim%2==0) fXdim +=1;
  if(fYdim%2==0) fYdim +=1;
  if(fZdim%2==0) fZdim +=1;
  
  ScalarField smoothFilter;
  
  smoothFilter.init(fXdim,fYdim,fZdim);
  smoothFilter.setVoxelSize(0,Xres); 
  smoothFilter.setVoxelSize(1,Yres); 
  smoothFilter.setVoxelSize(2,Zres);
 
  smoothFilter.ComputeSmoothingFilter(sigma);
    
  // assume kernel sizes are odd
  
  float *kernel_Data = smoothFilter.getAt(0,0,0);
  int order_x= (fXdim-1)/2;
  int order_y= (fYdim-1)/2;
  int order_z= (fZdim-1)/2;

  int ix_template, iy_template, iz_template;
  int ix1, iy1, iz1, ix2, iy2, iz2;

  // initialize to zero
  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  for(j=0; j<6; j++)
	    val[j] = 0.0f;
	}

  // create a unity mask for normalization
  
  ScalarField mask_float;
  mask_float.init(Xdim,Ydim,Zdim);
  mask_float.setVoxelSize(0,Xres); 
  mask_float.setVoxelSize(1,Yres); 
  mask_float.setVoxelSize(2,Zres);

  ScalarField mask_float_smooth;
  mask_float_smooth.init(Xdim,Ydim,Zdim);
  mask_float_smooth.setVoxelSize(0,Xres); 
  mask_float_smooth.setVoxelSize(1,Yres); 
  mask_float_smooth.setVoxelSize(2,Zres);


  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  float *mask_val = mask_float.getAt(ix,iy,iz);	  
	  mask_val[0] = (float) 1.0;
	}
	

  // (ix,iy,iz) is now a position in subject domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{

	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *warpField_val = warpField->getAt(ix,iy,iz);
	  ix_template = rint(warpField_val[0]);
	  iy_template = rint(warpField_val[1]);
	  iz_template = rint(warpField_val[2]);
	  
	  if(ix_template < 0 || ix_template >= Xdim || iy_template< 0 || iy_template >= Ydim || iz_template< 0 || iz_template>=Zdim)
	     continue;

	  smoothFilter.ComputeSmoothingFilter(sigma, warpField_val[0]-ix_template + order_x, warpField_val[1]-iy_template + order_y, warpField_val[2]-iz_template + order_z );    


	  float *mask_float_val= mask_float.getAt(ix,iy,iz);

	  // do trunc-Gaussian interpolation
	  for(ix1=0; ix1<fXdim; ix1++)
	    for(iy1=0; iy1<fYdim; iy1++)
	      for(iz1=0; iz1<fZdim; iz1++)
		{
		  int ind1 = iz1*fXdim*fYdim + iy1*fXdim + ix1;
		  
		  ix2 = ix_template -ix1 + order_x;
		  iy2 = iy_template -iy1 + order_y;
		  iz2 = iz_template -iz1 + order_z;
		  
		  if(ix2 < 0 || ix2 >= Xdim || iy2< 0 || iy2 >= Ydim || iz2< 0 || iz2>=Zdim)
		    continue;


		  float *val= dti_out->getAt(ix2,iy2,iz2);
		  for(j=0; j<6; j++)
		    val[j] += kernel_Data[ind1]*m_Data[ind*6+j];


		  float *mask_float_smooth_val= mask_float_smooth.getAt(ix2,iy2,iz2);

		  mask_float_smooth_val[0] += kernel_Data[ind1]*mask_float_val[0];

		}	  
	  
	}

  // point-wise normalization

  // (ix,iy,iz) is a position in template domain
  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  
	  float *mask_float_smooth_val= mask_float_smooth.getAt(ix,iy,iz);
	  
	  if(mask_float_smooth_val[0]>0)
	  for(j=0; j<6; j++)
	    val[j] /= mask_float_smooth_val[0];
	}
  
  // drop all tensors where trace falls below a threshold

  float *test_trace_val = dti_out->getAt(Zdim/2,Ydim/2,Zdim/2);
  float trace_cen = test_trace_val[0] + test_trace_val[1] + test_trace_val[2] ;

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  float *val= dti_out->getAt(ix,iy,iz);
	  
	  float trace = val[0]+val[1]+val[2];

	  if(trace < 0.01*trace_cen)
	    
	  for(j=0; j<6; j++)
	    val[j] = 0.0f;

	}


}

