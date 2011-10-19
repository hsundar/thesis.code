#include <iostream>
#include <QtCore>

#include "OptimizerOctSim.h"
#include "Volume.h"
#include "colors.h"


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << RED"Usage: "GRN << argv[0] << YLW" target_image moving_image"NRM << std::endl;
    std::cerr << std::endl << "Octree based Multi-Modality Registration" << std::endl \
              << "The Analyze and MetaIO formats are supported" << std::endl;

    std::cerr << "Hari Sundar, University of Pennsylvania, 28 February 2007." << std::endl;
    return -1;
  }

  // The filenames ...
  QString targetFilename(argv[1]);
  QString sourceFilename(argv[2]);

  // First read in the two images ...
  Volume *source = new Volume();
  Volume *target = new Volume();

  if (targetFilename.endsWith("mhd"))
    target->InitMhd(targetFilename);
  else if ( targetFilename.endsWith("hdr"))
    target->InitAnalyze(targetFilename);
  else {
    std::cerr << RED"Error: "NRM"Cannot load volume file "GRN << qPrintable(targetFilename) << NRM <<std::endl;
    return -2;
  }

  if (sourceFilename.endsWith("mhd"))
    source->InitMhd(sourceFilename);
  else if ( sourceFilename.endsWith("hdr"))
    source->InitAnalyze(sourceFilename);
  else {
    std::cerr << RED"Error: "NRM"Cannot load volume file "GRN << qPrintable(sourceFilename) << NRM <<std::endl;
    return -2;
  }

  // Now declare the oct_sim class ...
  oct_sim *sim = new oct_sim();
   
  sim->setHistBits(5);
  sim->initHists();

  sim->setTargetImage(target);
  sim->setSourceImage(source);

  // generate the Octree ...

  if (true) {
    sim->useOctree(true);
    target->SetMaxDepth(8);
    std::cout << "Computing Octree" << std::endl;
    target->ComputeLinearOctree();
    std::cout << "done Computing Octree" << std::endl;
  }

  // The optimizer ...
  const int dims = 6;

  double params[dims] = { 5, 5, 5, 0.1, 0.1, 0.1 };
  double steps[dims] = { 1, 1, 1, 0.01, 0.01, 0.01 };
  
  //OptimizerOctGD *optim = new OptimizerOctGD();
  OptimizerOctPB *optim = new OptimizerOctPB();
  optim->setInitialParameters (params, steps, dims);
  optim->setSimilarityMeasure(sim);
  // optim->init();

  std::cout << YLW"Starting Registration"NRM << std::endl;
  optim->optimize();
  std::cout << GRN"Finished Registration"NRM << std::endl;

  double *x = optim->getParameters ();
  std::cout << GRN"Final Params := ["YLW << x[0] <<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<", "<<x[4]<<", "<<x[5] \
            << GRN"] -> "BLU << optim->getCostFunctionValue () << NRM << std::endl;
  printf ("%i iterations, %i evaluations\n",
      optim->getFinalNumberOfIterations (),
      optim->getFinalNumberOfEvaluations ());

  delete optim;
  delete sim;
  delete source;
  delete target;

  return 0;
}

