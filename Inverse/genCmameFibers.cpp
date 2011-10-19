#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " Ns Nt out_File_Prefix forceFactor" << std::endl;
    return -1;
  }

  unsigned int sz = atoi(argv[1]);
  unsigned int nt = atoi(argv[2]);
  unsigned int nsz = sz;
  double ffac = 2.0;

  if (argc > 4)
    ffac = atof(argv[4]);

  std::cout << "Size is " << sz << std::endl;

  unsigned char *vol = new unsigned char [sz*sz*sz];

  for (unsigned int i=0; i<sz*sz*sz; i++)
    vol[i] = 0;


  char hdrname[256];	
  char volFile[256];
  sprintf(volFile, "%s.%d.img", argv[3], sz);

  //double *fldTemplate = new double [3*nsz*nsz*nsz];
  float *fldTemplate = new float [3*nsz*nsz*nsz];
  for (unsigned int i=0; i<3*nsz*nsz*nsz; i++)
    fldTemplate[i] = 0.0;

  unsigned int ht = 4*sz/5;
  unsigned int rm1 = sz/2.5;
  unsigned int rm2 = sz/4.5;
  int ctr = sz/2;
	 unsigned int ht2 = ht - rm1 + rm2;
	 
  unsigned int zmin, zmax;
  zmin = sz/10; zmax = zmin+ht;

  for (int k=0; k<ht; k++) {
    for (int j=0; j<sz; j++) {
      for (int i=0; i<sz; i++) {
        int _in=0;
        float fac = sqrt((1.0/ht)*k);
        unsigned int r1 = fac*rm1;
  
        double rr = sqrt((i-ctr)*(i-ctr) + (j-ctr)*(j-ctr));

        if ( rr <= r1 )
          _in = 1;
        
        if (k > (rm1 - rm2)) {
          fac  = sqrt((1.0/ht2)*(k-rm1+rm2));  
          unsigned int r2 = fac*rm2;

          if ( rr <= r2 )
            _in = 0;
        }

        if ( _in )
          vol [ sz*(j + (k+zmin)*sz) + i ] =  (unsigned char)200;
      } // i
    } // j
  } // k


  std::ofstream out(volFile);
  out.write((char *)vol, sz*sz*sz);
  out.close();

  sprintf(hdrname, "%s.%d.mhd", argv[3], sz);
  out.open(hdrname);
  out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
  out << "ElementSpacing = 2.816 2.816 2.816" <<  std::endl;
  out << "DimSize = " << sz << " " << sz << " " << sz << std::endl;
  out << "ElementType = MET_UCHAR" << std::endl;
  out << "ElementDataFile = " << volFile << std::endl;
  out.close();

  double *fld = new double [nsz*nsz*nsz];
  // Generate the fibers ...

  for (int k=0; k<ht; k++) {
    for (int j=0; j<nsz; j++) {
      for (int i=0; i<nsz; i++) {
        int _in=0;
        float fac = sqrt((1.0/ht)*k);
        unsigned int r1 = fac*rm1;
  
        double rr = sqrt((i-ctr)*(i-ctr) + (j-ctr)*(j-ctr));

        if ( rr <= r1 )
          _in = 1;
        
        if (k > (rm1 - rm2)) {
          fac  = sqrt((1.0/ht2)*(k-rm1+rm2));  
          unsigned int r2 = fac*rm2;

          if ( rr <= r2 )
            _in = 0;
        }

        if ( _in ) {
          fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)] = -(j-ctr); fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)] /= rr;
          fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+1] = i-ctr; fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+1] /= rr;
          double fac2 = 2*(rr-rm2)/(rm1-rm2);
          fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+2] = fac2 - 1.0;
          // Normalize
					double mag = sqrt( fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)]*fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)] +
														 fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+1]*fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+1] +
														 fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+2]*fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+2] );
					if (mag > 0.001) {
						fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)] /= mag;
						fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+1] /= mag;
						fldTemplate [ 3*(nsz*(j + (k+zmin)*nsz) + i)+2] /= mag;
					}
					fld[nsz*(j + (k+zmin)*nsz) + i] = 1000.0;
        }					
      } // i
    } // j
  } // k

  sprintf(volFile, "%s.%d.fibers", argv[3], sz);
  out.open(volFile);
  out.write((char *)fldTemplate, 3*nsz*nsz*nsz*sizeof(float));
  // out.write((char *)fldTemplate, 3*nsz*nsz*nsz*sizeof(double));
  out.close();

  sprintf(hdrname, "%s.%d.fibers.mhd", argv[3], sz);
  out.open(hdrname);
  out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
  out << "ElementSpacing = 2.816 2.816 2.816" <<  std::endl;
  //out << "ElementSpacing =  2.73067 2.73067 2.73067" << std::endl;
  out << "DimSize = " << nsz << " " << nsz << " " << nsz << std::endl;
  out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_FLOAT" << std::endl;
  out << "ElementDataFile = " << volFile << std::endl;
  out.close();

  // Now will write out the activation fields ...


  double *tau = new double[sz*sz*sz];

  char fldFile[256];
  unsigned int endSys = nt/3;

  for (unsigned int i=0; i<nt; i++) {
    double decay =0.0;
    if (i<endSys) 
      decay = ((double)i)/endSys;
    else
      decay = 1.0 - ((double)i - endSys)/(2*endSys);
    for (unsigned int j=0; j<sz*sz*sz; j++)
      tau[j] = decay*fld[j];
    sprintf(fldFile, "%s.%d.%.3d.fld", argv[3], sz, i);

    std::ofstream ofile(fldFile);
    ofile.write((char *)tau, nsz*nsz*nsz*sizeof(double));
    ofile.close();

    sprintf(hdrname, "%s_tau_%d.%.3d.mhd", argv[3], sz, i);
    out.open(hdrname);
    out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
    //out << "ElementSpacing =  2.73067 2.73067 2.73067" << std::endl;
    out << "ElementSpacing = 2.816 2.816 2.816" <<  std::endl;
    out << "DimSize = " << nsz << " " << nsz << " " << nsz << std::endl;
    out << "ElementType = MET_DOUBLE" << std::endl;
    out << "ElementDataFile = " << fldFile << std::endl;
    out.close();
  }


  return 0;
}
