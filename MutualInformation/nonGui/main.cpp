#include <iostream> 
#include "oct_sim.h"
#include "TransMatrix.h"


int main(int argc, char **argv) {

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " sourceImg targetImg useOct" << std::endl;
    return -1;
  } // if (argc < 3)

  bool useOct = false;
  if (argc>3)
    useOct = true;

  QString srcName(argv[1]);
  QString tarName(argv[2]);


  Volume *source = new Volume();
  Volume *target = new Volume();

  if ( srcName.endsWith("hdr") ) {
    source->InitAnalyze(srcName);
  } else if ( srcName.endsWith("mhd") ) {
    source->InitMhd(srcName);
  } else {
    std::cerr << "Unknown file format for source image" << std::endl;
    return -1;
  }

  if ( tarName.endsWith("hdr") )
    target->InitAnalyze(tarName);
  else if ( tarName.endsWith("mhd") )
    target->InitMhd(tarName);
  else {
    std::cerr << "Unknown file format for target image" << std::endl;
    return -1;
  }

  // set up the similarity evaluator ...
  oct_sim *sim = new oct_sim();

  sim->setSourceImage(source);
  sim->setTargetImage(target);
  sim->setHistBits(6);
  sim->initHists();


  if (useOct) {
    target->SetMaxDepth(8);
    target->ComputeLinearOctree();
    // std::cout << "Size of Octree is " << target->octSize() << std::endl;
    sim->useOctree(true);
  }

  TransMatrix mat, mat1, mat2; 

  int iter =10;
  std::cout << "STARTING " << std::endl << std::endl;
  for (int i =-iter; i < iter; i++) { 
    std::cout << "Iteration " << i; // << std::endl;
    mat = TransMatrix::Translation(5.*i, 0., 0. );
    // mat = TransMatrix::Translation(0., 5.*i, 0. );
    // mat = TransMatrix::Translation(0., 0., 2.*i );
    //mat1 = TransMatrix::RotationX(0.1*i);
    //mat2 = TransMatrix::Translation(-128, -128,-128);
    //mat1 = mat1.MultMatrixLeftBy(mat2);
    //mat2 = TransMatrix::Translation(128, 128,128);
    // mat2 = TransMatrix::Translation(4*i, 2*i, 0.5*i);
    //mat = mat1.MultMatrixRightBy(mat2);
    double ent = sim->getNMI(mat.GetDataPointer());
    std::cout << " " <<  ent << std::endl;
  } 

  delete sim;
  delete source;
  delete target;

  return 0;
}

