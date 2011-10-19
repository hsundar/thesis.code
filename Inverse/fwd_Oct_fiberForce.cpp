static char help[] = "Driver for a linear elastodynamic problem (Hyperbolic)";

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "petscksp.h"
#include "petscda.h"

#include "parUtils.h"

#include "timeInfo.h"
#include "feMatrix.h"
#include "feVector.h"
#include "femUtils.h"
#include "timeStepper.h"
#include "newmark.h"
#include "elasStiffness.h"
#include "elasMass.h"
#include "raleighDamping.h"
#include "cardiacForce.h"

float uniform() {
  return float(rand()) / RAND_MAX; // [0,1)
}

double gaussian(double mean = 0.5, double std_deviation = 0.1) {
  static double t = 0;
  double x1, x2, r;

  // reuse previous calculations
  if (t) {
    const double tmp = t;
    t = 0;
    return mean + std_deviation * tmp;
  }

  // pick randomly a point inside the unit disk
  do {
    x1 = 2 * uniform() - 1;
    x2 = 2 * uniform() - 1;
    r = x1 * x1 + x2 * x2;
  } while (r >= 1);

  // Box-Muller transform
  r = sqrt(-2.0 * log(r) / r);

  // save for next call
  t = r * x2;

  // only use one of the coordinates of a bivariate distribution
  return mean + std_deviation * r * x1;
}

/*
double randgauss( double min, double max, double sigma, double centre) {
  double random = (min + (max-min) * (double)rand()/RAND_MAX); //create random domain between [min,max]

  double tmp = (random-centre)/sigma;
  double gauss= exp(-tmp*tmp/2); //gaussian formula

  return gauss;
}
*/

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, "cmame.opt", help);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // The global domain size
  double gSize[3];
  gSize[0] = 1.; gSize[1] = 1.; gSize[2] = 1.;

  // Parameters for the balancing algorithm.
  // Refer to manual for details ... 
  bool incCorner = 1; // balance across corners = true  
  unsigned int maxNumPts= 1; // maximum number of points per octant
  unsigned int dim=3; // spatial dimensions 

  // unsigned int maxDepth=8; // maximum depth of the octree, has to be <= 30
  int maxDepth=8; // maximum depth of the octree, has to be <= 30

  int Ns = 32;
  unsigned int dof = 3;

  char problemName[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];

  double t0 = 0.0;
  double dt = 0.1;
  double t1 = 1.0;

  Vec rho;        // density - elemental scalar
  Vec lambda;     // Lame parameter - lambda - elemental scalar
  Vec mu;         // Lame parameter - mu - elemental scalar
  Vec fibers;     // Fiber orientations - nodal vector (3-dof)

  std::vector<Vec> tau;        // the scalar activation - nodal scalar

  std::vector<ot::TreeNode> linOct, balOct, newLinOct;
  std::vector<double> pts;

  // Initial conditions
  Vec initialDisplacement; 
  Vec initialVelocity;

  double nu, E;
  nu = 0.45;
  E = 1000;

  timeInfo ti;

  PetscTruth mf = PETSC_FALSE;
  bool mfree = false;

  PetscOptionsGetTruth(0, "-mfree", &mf, 0);

  if (mf == PETSC_TRUE) {
    mfree = true;
  } else
    mfree = false;

  double ctrst = 1.0;
  // get Ns
  CHKERRQ ( PetscOptionsGetInt(0,"-Ns",&Ns,0) );
  CHKERRQ ( PetscOptionsGetInt(0,"-mdepth",&maxDepth,0) );

  CHKERRQ ( PetscOptionsGetScalar(0,"-ctrst",&ctrst,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t0",&t0,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-nu",&nu,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-Youngs",&E,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t1",&t1,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-dt",&dt,0) );
  CHKERRQ ( PetscOptionsGetString(PETSC_NULL,"-pn",problemName,PETSC_MAX_PATH_LEN-1,PETSC_NULL));

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  if (!rank) {
    std::cout << "Grid size is " << Ns+1 << " and NT is " << (int)ceil(1.0/dt) << std::endl;
    std::cout << "MaxDepth is " << maxDepth << std::endl;
  }

  /*
  for (int i=0; i<3*Ns*Ns*Ns; i++) {
    double val = gaussian(); //randgauss(0., 1., 2.0, 0.);
    // std::cout << val << std::endl;
    pts.push_back(val);
  }

  MPI_Barrier(MPI_COMM_WORLD);  
*/
  /*********************************************************************** */
  // CONSTRUCT: Construct the linear octree from the points ...
  /*********************************************************************** */
  // ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI_COMM_WORLD);

  // The points are not needed anymore, and can be cleared to free memory.
  // pts.clear();

  if (!rank) {
    ot::readNodesFromFile("test.256.oct", newLinOct);
    std::cout << "Finished reading" << std::endl;
  }

  // std::sort(linOct.begin(), linOct.end());
  std::cout << rank << " Original octree size is " << newLinOct.size() << std::endl;

  /*
  par::Partition<ot::TreeNode>(linoct, newLinOct, MPI_COMM_WORLD);
  linOct.clear();
  */
  par::sampleSort<ot::TreeNode>(newLinOct, linOct, MPI_COMM_WORLD);
  newLinOct.clear();

  std::cout << rank << ": after Part octree size is " << linOct.size() << std::endl;

  /*********************************************************************** */
  // BALANCE: Balance the linear octree to enforce the 2:1 balance conditions.
  /*********************************************************************** */
  ot::balanceOctree (linOct, balOct, dim, maxDepth, incCorner, MPI_COMM_WORLD);

  std::cout << "Balanced octree size is " << balOct.size() << std::endl;

  // The linear octree (unbalanced) can be cleared to free memory.
  linOct.clear();

  // If desired the octree can be written to a file using the supplied routine ...
  ot::writeNodesToFile("filename.oct", balOct);

  /*********************************************************************** */
  // MESH : Construct the octree-based Distruted Array.
  /*********************************************************************** */
  ot::DA da(balOct,MPI_COMM_WORLD);
  balOct.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  if (!rank)
    std::cout <<"Finshed Meshing" << std::endl;

  PetscFinalize();
  return 0;

  // create Matrices and Vectors
  elasMass *Mass = new elasMass(feMat::OCT); // Mass Matrix
  elasStiffness *Stiffness = new elasStiffness(feMat::OCT); // Stiffness matrix
  raleighDamping *Damping = new raleighDamping(feMat::OCT); // Damping Matrix

  cardiacForce *Force = new cardiacForce(feVec::OCT); // Force Vector

  // create vectors 

  da.createVector(rho, false, false, 1);
  da.createVector(mu, false, false, 1);
  da.createVector(lambda, false, false, 1);

  da.createVector(initialDisplacement, false, false, dof);
  da.createVector(initialVelocity, false, false, dof);

  // Set initial conditions
  CHKERRQ( VecSet ( initialDisplacement, 0.0) ); 
  CHKERRQ( VecSet ( initialVelocity, 0.0) );

  // Homogeneous material properties ...
  /*
  CHKERRQ( VecSet ( rho, 1.0) ); 

  nu = 0.45; E = 1000;
  double mmu = E/(2*(1+nu));
  double llam = E*nu/((1+nu)*(1-2*nu));

  CHKERRQ( VecSet ( mu, mmu) ); 
  CHKERRQ( VecSet ( lambda, llam) );
  */


  PetscScalar *muArray, *lamArray, *rhoArray;
  // Read in material properties from file ...
  unsigned int elemSize = Ns*Ns*Ns;

  unsigned char *tmp_mat = new unsigned char[elemSize];  
  double *tmp_tau = new double[elemSize];
  double *tmp_fib = new double[dof*elemSize];

  // generate filenames & read in the raw arrays first ...
  std::ifstream fin;

  sprintf(filename, "%s.%d.img", problemName, Ns); 
  fin.open(filename, std::ios::binary); fin.read((char *)tmp_mat, elemSize); fin.close();

  nu = 0.35; E = 1000;
  double mmu = E/(2*(1+nu));
  double llam = E*nu/((1+nu)*(1-2*nu));
  nu = 0.45; E = 1000*ctrst;
  double mmu2 = E/(2*(1+nu));
  double llam2 = E*nu/((1+nu)*(1-2*nu));

  da.vecGetBuffer(mu, muArray,true,true,false,1);
  da.vecGetBuffer(lambda, lamArray,true,true,false,1);
  da.vecGetBuffer(rho, rhoArray,true,true,false,1);

  for ( da.init<ot::DA::ALL>(), da.init<ot::DA::WRITABLE>(); da.curr() < da.end<ot::DA::ALL>(); da.next<ot::DA::ALL>()) {
    unsigned int i = da.curr();
    Point pt;
    pt = da.getCurrentOffset();

    int indx = pt.z()*Ns*Ns + pt.y()*Ns + pt.x();

    if ( tmp_mat[indx] ) {
      muArray[i] = mmu2;
      lamArray[i] = llam2;
      rhoArray[i] = 1.0;
    } else {
      muArray[i] = mmu;
      lamArray[i] = llam;
      rhoArray[i] = 1.0;
    }
  }

  da.vecRestoreBuffer(mu, muArray,true,true,false,1);
  da.vecRestoreBuffer(lambda, lamArray,true,true,false,1);
  da.vecRestoreBuffer(rho, rhoArray,true,true,false,1);

  delete [] tmp_mat;

  // Now set the activation ...
  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));

  Vec tauVec, tmpTau;

  // Create elemental vector ...
  da.createVector(tmpTau, true, true, 1);
  da.createVector(fibers, true, true, 3);
  
  PetscScalar *tauArray;

  // load the fibers ...
  da.vecGetBuffer(fibers, tauArray, true, true, 3);
    
  sprintf(filename, "%s.%d.fibers", problemName, Ns);
  std::ifstream fin3(filename, std::ios::binary); fin3.read((char *)tmp_fib, dof*elemSize*sizeof(double)); fin3.close();

  for ( da.init<ot::DA::ALL>(), da.init<ot::DA::WRITABLE>(); da.curr() < da.end<ot::DA::ALL>(); da.next<ot::DA::ALL>()) {
    unsigned int i = da.curr();
    Point pt;
    pt = da.getCurrentOffset();

    int indx = pt.z()*Ns*Ns + pt.y()*Ns + pt.x();

    for (int d=0; d<dof; d++) {
      tauArray[dof*i+d] = tmp_tau[dof*indx+d];
    }
  }

  da.vecRestoreBuffer(fibers, tauArray, true, true, 3);

  delete [] tmp_fib;


  // loop through time steps
  for (unsigned int t=0; t<numSteps+1; t++) {
    // a. Create new nodal vector ...
    da.createVector(tauVec, false, true, 1);

    VecSet( tmpTau, 0.0 );
    da.vecGetBuffer(tmpTau, tauArray, true, true, 1);

    // b. read in the activation 
    // std::cout << "Setting force vectors" << std::endl;
    sprintf(filename, "%s.%d.%.3d.fld", problemName, Ns, t);
    // std::cout << "Reading force file " << filename << std::endl;
    fin.open(filename); fin.read((char *)tmp_tau, elemSize*sizeof(double)); fin.close();

    // c. set the values ...
    for ( da.init<ot::DA::ALL>(), da.init<ot::DA::WRITABLE>(); da.curr() < da.end<ot::DA::ALL>(); da.next<ot::DA::ALL>()) {
      unsigned int i = da.curr();
      Point pt;
      pt = da.getCurrentOffset();

      int indx = pt.z()*Ns*Ns + pt.y()*Ns + pt.x();
      
      tauArray[i] = tmp_tau[indx];
    }

    // restore
    da.vecRestoreBuffer(tmpTau, tauArray, true, true, 1);
    // d. element2node
    elementToNode(da, tmpTau, tauVec, 1);

    // store in vector 
    tau.push_back(tauVec);
  }

  // Setup Matrices and Force Vector ...

  Mass->setProblemDimensions(1.0, 1.0, 1.0);
  Mass->setDA(&da);
  Mass->setDof(dof);
  Mass->setDensity(rho);

  Stiffness->setProblemDimensions(1.0, 1.0, 1.0);
  Stiffness->setDA(&da);
  Stiffness->setDof(dof);
  Stiffness->setLame(lambda, mu);

  Damping->setAlpha(0.0);
  Damping->setBeta(0.00075);
  Damping->setMassMatrix(Mass);
  Damping->setStiffnessMatrix(Stiffness);
  Damping->setDA(&da);
  Damping->setDof(dof);

  // Force Vector
  Force->setProblemDimensions(1.0,1.0,1.0);
  Force->setDA(&da);
  // Force->setFDynamic(tau);
  Force->setActivationVec(tau);
  Force->setFiberOrientations(fibers);
  Force->setTimeInfo(&ti);

  // Newmark time stepper ...
  newmark *ts = new newmark; 

  ts->setMassMatrix(Mass);
  ts->setDampingMatrix(Damping);
  ts->setStiffnessMatrix(Stiffness);
  ts->damp(false);
  ts->setTimeFrames(1);
  ts->storeVec(false);

  ts->setForceVector(Force);

  ts->setInitialDisplacement(initialDisplacement);
  ts->setInitialVelocity(initialVelocity);

  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward
  ts->useMatrixFree(mfree);

  //if (!rank)
  //  std::cout << RED"Initializing Newmark"NRM << std::endl;
  double itime = MPI_Wtime();
  ts->init(); // initialize IMPORTANT 
  //if (!rank)
  //  std::cout << RED"Starting Newmark Solve"NRM << std::endl;
  double stime = MPI_Wtime();
  ts->solve();// solve 
  double etime = MPI_Wtime();
  //if (!rank)
  //  std::cout << GRN"Done Newmark"NRM << std::endl;
  if (!rank) {
    std::cout << "Total time for init is " << stime - itime << std::endl;
    std::cout << "Total time for solve is " << etime - stime << std::endl;
  }

  PetscFinalize();
}

