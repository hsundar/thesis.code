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
  double ffac = 2.0;

  if (argc > 4)
    ffac = atof(argv[4]);

  std::cout << "Size is " << sz << std::endl;

  unsigned char *vol = new unsigned char [sz*sz*sz];

  for (unsigned int i=0; i<sz*sz*sz; i++)
    vol[i] = 0;


  char volFile[256];
  sprintf(volFile, "%s.%d.img", argv[3], sz);

  double *fldTemplate = new double [3*sz*sz*sz];
  for (unsigned int i=0; i<3*sz*sz*sz; i++)
    fldTemplate[i] = 0.0;

  unsigned int ht = 3*sz/5;
  unsigned int rm1 = sz/3;
  unsigned int rm2 = sz/4;
  int ctr = sz/2;
  unsigned int ht2 = ht - rm1 + rm2;


  unsigned int zmin, zmax;
  zmin = sz/5; zmax = zmin+ht;

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

        if ( !_in ) // outside ...
          vol [ sz*(j + (k+zmin)*sz) + i ] =  (unsigned char)0;
        else {
          vol [ sz*(j + (k+zmin)*sz) + i ] =  (unsigned char)200;
          if (rr > 0.00001) {
            fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) ] = (-(double)(ctr-i))/rr;
            fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) + 1] = (-(double)(ctr-j))/rr;
            /*
	    if ( isnan(fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) ]) ) {
              std::cout << "Have NAN" << std::endl;
              std::cout << "RR is " << rr << std::endl;
            } */

          } else {
            fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) ] = 0.0;
            fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) + 1] = 0.0;
          }
          fldTemplate [ 3*( sz*(j + (k+zmin)*sz) + i ) + 2] = 0.0;
        }
      } // i
    } // j
  } // k

  std::ofstream out(volFile);
  out.write((char *)vol, sz*sz*sz);
  // close the file ...
  out.close();


  // Now will write out the force fields ...
  double *fld = new double [3*sz*sz*sz];
  
  for (unsigned int i=0; i<3*sz*sz*sz; i++)
    fld[i] = 0.0;

  char fldFile[256];
  unsigned int endSys = nt/3;

  for (unsigned int i=0; i<nt; i++) {
    sprintf(fldFile, "%s.%d.%.3d.fld", argv[3], sz, i);
    
    double decay=1.0;
    if ( i <= endSys ) {
      decay = ((double)i)/endSys;
    } else {
      decay = ((double)(nt - i))/(2.0*endSys);
    }

    for (unsigned int i=0; i<3*sz*sz*sz; i++)
      fld[i] = decay*ffac*fldTemplate[i];

    std::ofstream ofile(fldFile);
    ofile.write((char *)fld, 3*sz*sz*sz*sizeof(double));
    ofile.close();
  }


  return 0;
}
