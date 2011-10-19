#include <fstream>
#include <iostream>
#include <algorithm>

#include "Point.h"
#include "Field.h"
#include "colors.h"
/*
   Takes the fiber orientation file as input and generates the following
   - image file at desired resolution
   - fiber orientation at desired resolution
   - activation field at the desired resolution

   Assumes that the input fiber orientation will be provided in the MetaIO format.

*/	

int main(int argc, char** argv) {
  if(argc < 4) {
    std::cerr << RED"Usage: "NRM << argv[0] << " fiber_orientation[MHD] grid_size[uint] timesteps[uint] problem_name[string]" << std::endl;
    std::cerr << std::endl;
    std::cerr << GRN"\texample: "NRM << argv[0] << " fibers.mhd 32 100 heart" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Takes the fiber orientation file as input and generates the following" << std::endl;
    std::cerr << "\t - image file at desired resolution" << std::endl;
    std::cerr << "\t - fiber orientation at desired resolution" << std::endl;
    std::cerr << "\t - activation fields at the desired resolution" << std::endl;
    std::cerr << std::endl << "Assumes that the input fiber orientation will be provided in the MetaIO format." << std::endl;

    return -1;
  }

  Field<float, 3> fibers;
  Field<float, 1> fa;
  unsigned int Ns;
  unsigned int Nt;

  char fname[256];
  char hdrname[256];

  // The inputs ...
  fibers.init(argv[1]);
  // fa.init("c09.mhd");
  Ns = atoi(argv[2]);
  Nt = atoi(argv[3]);

  std::cout << "Read in Fibers of size ("GRN << fibers.getSize(0) << NRM","GRN << fibers.getSize(1) << NRM","GRN << fibers.getSize(2) << NRM")";
  std::cout << " and spacing ("GRN << fibers.getVoxelSize(0) << NRM","GRN << fibers.getVoxelSize(1) << NRM","GRN << fibers.getVoxelSize(2) << NRM")" << std::endl;

  /*  We shall generate the following ... in the order
      1. Subsampled and isotropic fiber orientations ...
      2. Image file from the new fibers
      3. activation for the image ...
      */

  int x,y,z;
  x = fibers.getSize(0);
  y = fibers.getSize(1);
  z = fibers.getSize(2);
  double sx, sy, sz;
  sx = fibers.getVoxelSize(0);
  sy = fibers.getVoxelSize(1);
  sz = fibers.getVoxelSize(2);

  double lx, ly, lz, maxl;
  lx = sx*x;
  ly = sy*y;
  lz = sz*z;

  maxl = std::max(lx, lz); // assuming that lx=ly

  // Allocate memory ...
  // double *fiber_new  = new double[3*Ns*Ns*Ns];
  float *fiber_new  = new float[3*Ns*Ns*Ns];
  // float *fa_new  = new float[Ns*Ns*Ns];

  unsigned char *geo = new unsigned char[Ns*Ns*Ns];

  /***** First the new fibers *****/

  // Determine the padding required to make it isotropic and equi-dimensional
  double padZ = 0.; //.5*(maxl - lz);
  // assert (padZ >= 0.0); // for now assuming Z is short 
  for (unsigned int k=0; k<Ns; k++) {
    for (unsigned int j=0; j<Ns; j++) {
      for (unsigned int i=0; i<Ns; i++) {
        Point pt;
        // compute the correct index to sample from.
        double px, py, pz;
        px = (maxl*i)/Ns/sx;
        py = (maxl*j)/Ns/sy;
        //pz = ((maxl*k)/Ns - padZ)/sz;
        pz = (maxl*k)/Ns/sz;

        pt = fibers.getAt(px, py, pz);
        pt.normalize();
        // set the values ...
        fiber_new[3*((k*Ns+j)*Ns+i)] 		= pt.x();
        fiber_new[3*((k*Ns+j)*Ns+i)+1] 	= pt.y();
        fiber_new[3*((k*Ns+j)*Ns+i)+2] 	= pt.z();

        // float * val = fa.getAt((int)px, (int)py, (int)pz);
        // fa_new[(k*Ns+j)*Ns+i] = *val;

        if (pt.abs() > 0.01)
          geo[(k*Ns+j)*Ns+i] = 200;
        else
          geo[(k*Ns+j)*Ns+i] = 0;
      }
    }
  }		


  double *tau = new double[Ns*Ns*Ns];
  unsigned char *uctau = new unsigned char[Ns*Ns*Ns];

  std::cout << "Done generating fibers" << std::endl;
  /*
  // generate the nodal fibers ...
  delete [] fiber_new;
  fiber_new = new double[3*(Ns+1)*(Ns+1)*(Ns+1)];
  double *tau = new double[(Ns+1)*(Ns+1)*(Ns+1)];
  // unsigned char *uctau = new unsigned char[(Ns+1)*(Ns+1)*(Ns+1)];

  for (unsigned int k=0; k<Ns+1; k++) {
  for (unsigned int j=0; j<Ns+1; j++) {
  for (unsigned int i=0; i<Ns+1; i++) {
  Point pt;
  // compute the correct index to sample from.
  double px, py, pz;
  px = (maxl*i)/(Ns+1)/sx; // - maxl/Ns/sx/2;
  py = (maxl*j)/(Ns+1)/sy ; //- maxl/Ns/sx/2;
  pz = ((maxl*k)/(Ns+1) - padZ)/sz; // - maxl/Ns/sx/2;

  pt = fibers.getAt(px, py, pz);
  pt.normalize();
  // set the values ...
  fiber_new[3*((k*(Ns+1)+j)*(Ns+1)+i)] 		= pt.x();
  fiber_new[3*((k*(Ns+1)+j)*(Ns+1)+i)+1] 	= pt.y();
  fiber_new[3*((k*(Ns+1)+j)*(Ns+1)+i)+2] 	= pt.z();
  }
  }
  }
  */
  std::ofstream out;

  // Write out the activation ...
  double t1 = 0.3333*Nt;
  double fac= 200.0;
  for (unsigned t=0; t<Nt; t++) {
    for (unsigned int k=0; k<Ns; k++) {
      // determine the activation value for this slice ...
      double act=0.0;
      double sliceOnTime = (t1*(double)k)/(Ns);
      double sliceMidTime = sliceOnTime + t1; //(t2*(double)k)/(Ns+1);
      double sliceOffTime = sliceOnTime + 2*t1;
      if ( (t < sliceOnTime) || (t > sliceOffTime) )
        act =0.0;
      else if ( t < sliceMidTime) {
        double sliceTimeFac = (double)t - sliceOnTime;
        sliceTimeFac /= t1;
        act = fac*(sliceTimeFac*sliceTimeFac);
      } else {
        double sliceTimeFac = (double)t - sliceMidTime;
        sliceTimeFac = 1.0 - sliceTimeFac/t1;
        act = fac*(sliceTimeFac*sliceTimeFac);
      }

      for (unsigned int ij=0; ij<(Ns)*(Ns); ij++) {
        unsigned int idx = k*(Ns)*(Ns) + ij;
        double fiber_norm = sqrt(fiber_new[3*idx]*fiber_new[3*idx] + fiber_new[3*idx+1]*fiber_new[3*idx+1] + fiber_new[3*idx+2]*fiber_new[3*idx+2]);
        tau[idx] = act*fiber_norm;

        if ( (tau[idx]>0) && (tau[idx]<255) )
          uctau[idx] = (unsigned char)tau[idx];
        else 
          uctau[idx] = 0;
      }
    }
    // write out the file
    sprintf(fname, "%s.%d.%.3d.fld", argv[4], Ns, t);
    out.open(fname, std::ios::binary);
    out.write((char *)tau, (Ns)*(Ns)*(Ns)*sizeof(double));
    // out.write((char *)uctau, (Ns)*(Ns)*(Ns));
    out.close();
    sprintf(hdrname, "%s_tau_%d.%.3d.mhd", argv[4], Ns, t);
    out.open(hdrname);
    out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
    out << "ElementSpacing = " << maxl/(Ns) << " " << maxl/(Ns) << " " << maxl/(Ns) << std::endl;
    out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
    out << "ElementType = MET_DOUBLE" << std::endl;
    // out << "ElementType = MET_UCHAR" << std::endl;
    out << "ElementDataFile = " << fname << std::endl;
    out.close();

  }

  // write the files ...

  // FIBERS
  sprintf(fname, "%s.%d.fibers", argv[4], Ns);
  out.open(fname);
  // out.write((char *)fiber_new, (Ns)*(Ns)*(Ns)*3*sizeof(double));
  out.write((char *)fiber_new, (Ns)*(Ns)*(Ns)*3*sizeof(float));
  out.close();
  sprintf(hdrname, "%s.%d.fibers.mhd", argv[4], Ns);
  out.open(hdrname);
  out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
  out << "ElementSpacing = " << maxl/(Ns) << " " << maxl/(Ns) << " " << maxl/(Ns) << std::endl;
  out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
  // out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_DOUBLE" << std::endl;
  out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_FLOAT" << std::endl;
  out << "ElementDataFile = " << fname << std::endl;
  out.close();

  // FA
  /*
  sprintf(fname, "%s.%d.fa", argv[4], Ns);
  out.open(fname);
  out.write((char *)fa_new, (Ns)*(Ns)*(Ns)*sizeof(float));
  out.close();
  sprintf(hdrname, "%s.%d.fa.mhd", argv[4], Ns);
  out.open(hdrname);
  out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
  out << "ElementSpacing = " << maxl/(Ns) << " " << maxl/(Ns) << " " << maxl/(Ns) << std::endl;
  out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
  // out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_DOUBLE" << std::endl;
  out << "ElementNumberOfChannels = 1" << std::endl << "ElementType = MET_FLOAT" << std::endl;
  out << "ElementDataFile = " << fname << std::endl;
  out.close();
*/
  // IMAGE
  sprintf(fname, "%s.%d.img", argv[4], Ns);
  out.open(fname);
  out.write((char *)geo, Ns*Ns*Ns);
  out.close();
  sprintf(hdrname, "%s.%d.mhd", argv[4], Ns);
  out.open(hdrname);
  out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
  out << "ElementSpacing = " << maxl/Ns << " " << maxl/Ns << " " << maxl/Ns << std::endl;
  out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
  out << "ElementType = MET_UCHAR" << std::endl;
  out << "ElementDataFile = " << fname << std::endl;
  out.close();

  // clean up
  delete [] fiber_new;
  delete [] geo;
  delete [] tau;
  // delete [] uctau;
  return 0;
}
