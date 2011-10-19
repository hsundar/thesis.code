#include <iostream>
#include <fstream>

using namespace std;

int main (int argc, char ** argv) {
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " in_def_fld out_def_fld x y z" << endl;
		//                                     1           2     3 4 5
		return 1;
	}


	int x,y,z;
	x = atoi(argv[3]); y = atoi(argv[4]); z = atoi(argv[5]);
	
	float * A = new float[x*y*z*3];
	float * B = new float[x*y*z*3];

	// read in ....
	ifstream A_in(argv[1]);
	A_in.read((char *)A, x*y*z*3*sizeof(float));
	A_in.close();

	for (int k=0; k<z; k++) {
	  //	cout << "Processing slice " << k << endl;
		for (int j=0; j<y; j++)
		for (int i=0; i<x; i++) { 
			int ii = 3*(k*x*y + j*x + i);
			B[ii] = i + A[ii+1];
			B[ii+1] = j + A[ii];
			B[ii+2] = k + A[ii+2];

			if(i==x/2 && j==y/2 && k == z/2)
			  {
			    std::cout<<"\nCentral slice Out "<<B[ii]<<" "<<B[ii+1]<<" "<<B[ii+2]<<"\n";
			    std::cout<<"\nCentral slice In "<<A[ii]<<" "<<A[ii+1]<<" "<<A[ii+2]<<"\n";
			  }
		}
	}

	// write out ...
	ofstream out(argv[2]);
	out.write((char *)B, x*y*z*3*sizeof(float));
	out.close();

	delete [] A;
	delete [] B;

	return 0;
}

