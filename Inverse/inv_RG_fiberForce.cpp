static char help[] = "Driver for a estimate cardiac activations";

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "externVars.h"

#include "petscksp.h"
#include "petscda.h"

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
#include "parametricActivationInverse.h"
#include "radialBasis.h"
#include "bSplineBasis.h"

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

int main(int argc, char **argv)
{       
  PetscInitialize(&argc, &argv, "elas.opt", help);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int Ns = 32;
  unsigned int dof = 3;

  char problemName[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];

  double t0 = 0.0;
  double dt = 0.1;
  double t1 = 1.0;
  double beta = 0.000001;

  // double dtratio = 1.0;
  DA  da;         // Underlying scalar DA - for scalar properties
  DA  da3d;       // Underlying vector DA - for vector properties

  Vec rho;        // density - elemental scalar
  Vec lambda;     // Lame parameter - lambda - elemental scalar
  Vec mu;         // Lame parameter - mu - elemental scalar
  Vec fibers;     // Fiber orientations - elemental vector (3-dof)

  std::vector<Vec> tau;        // the scalar activation - nodal scalar

  // Initial conditions
  Vec initialDisplacement; 
  Vec initialVelocity;

  timeInfo ti;

  double nu, E;
  nu = 0.45;
  E = 1000;

  PetscTruth mf = PETSC_FALSE;
  bool mfree = false;

  PetscOptionsGetTruth(0, "-mfree", &mf, 0);

  if (mf == PETSC_TRUE) {
    mfree = true;
  } else 
    mfree = false;

  double ctrst = 10.0;

  int parFac = 2;
  int numParams = 120;

  CHKERRQ ( PetscOptionsGetInt(0,"-pFac", &parFac,0) );

  numParams = parFac*parFac*parFac*5;
  if (!rank)
    std::cout << "Total number of unknowns is " << numParams << std::endl;

  // get Ns
  CHKERRQ ( PetscOptionsGetInt(0,"-Ns",&Ns,0) );

  CHKERRQ ( PetscOptionsGetScalar(0,"-ctrst",&ctrst,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t0",&t0,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-nu",&nu,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-Youngs",&E,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t1",&t1,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-dt",&dt,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-beta",&beta,0) );
  CHKERRQ ( PetscOptionsGetString(PETSC_NULL,"-pn",problemName,PETSC_MAX_PATH_LEN-1,PETSC_NULL));

  if (!rank) {
    std::cout << "Problem size is " << Ns+1 << " spatially and NT = " << (int)ceil(1.0/dt) << std::endl;
  }

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  // create DA
  CHKERRQ ( DACreate3d ( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, 
        Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
        1, 1, 0, 0, 0, &da) );
  CHKERRQ ( DACreate3d ( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, 
        Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
        dof, 1, 0, 0, 0, &da3d) );

  elasMass *Mass = new elasMass(feMat::PETSC); // Mass Matrix
  elasStiffness *Stiffness = new elasStiffness(feMat::PETSC); // Stiffness matrix
  raleighDamping *Damping = new raleighDamping(feMat::PETSC); // Damping Matrix

  cardiacForce *Force = new cardiacForce(feVec::PETSC); // Force Vector

  // create vectors 
  CHKERRQ( DACreateGlobalVector(da, &rho) );
  CHKERRQ( DACreateGlobalVector(da, &mu) );
  CHKERRQ( DACreateGlobalVector(da, &lambda) );

  CHKERRQ( DACreateGlobalVector(da3d, &initialDisplacement) );
  CHKERRQ( DACreateGlobalVector(da3d, &initialVelocity) );

  // Set initial conditions
  CHKERRQ( VecSet ( initialDisplacement, 0.0) ); 
  CHKERRQ( VecSet ( initialVelocity, 0.0) );

  VecZeroEntries( mu );
  VecZeroEntries( lambda );
  VecZeroEntries( rho );

  int x, y, z, m, n, p;
  int mx,my,mz, xne, yne, zne;

  CHKERRQ( DAGetCorners(da, &x, &y, &z, &m, &n, &p) ); 
  CHKERRQ( DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ); 

  if (x+m == mx) {
    xne=m-1;
  } else {
    xne=m;
  }
  if (y+n == my) {
    yne=n-1;
  } else {
    yne=n;
  }
  if (z+p == mz) {
    zne=p-1;
  } else {
    zne=p;
  }

  double acx,acy,acz;
  double hx = 1.0/((double)Ns);

  // Generate the basis ...
  std::vector < radialBasis > spatialBasis;
  bSplineBasis temporalBasis(3, 5); // this creates and sets up the basis ...
  // temporalBasis.knot();

  double fac = 1.0/parFac;
  double ssq = fac*fac/5.5452; // 8log(2) = 5.5452
  PetscPrintf(0, "SSQ is %f\n", ssq); 
  // Now to set up the radial bases ...
  for (int k=0; k<parFac; k++) {
    for (int j=0; j<parFac; j++) {
      for (int i=0; i<parFac; i++) {
        // std::cout << "Adding radial basis at: " << fac/2+i*fac << ", " << fac/2+j*fac << ", " << fac/2+k*fac << std::endl;
        radialBasis tmp(Point( fac/2+i*fac,fac/2+j*fac,fac/2+k*fac), Point(ssq,ssq,ssq));
        spatialBasis.push_back(tmp);
      }
    }
  }

  // SET MATERIAL PROPERTIES ...

  // @todo - Write routines to read/write in Parallel
  // allocate for temporary buffers ...
  unsigned int elemSize = Ns*Ns*Ns;
  unsigned int nodeSize = (Ns+1)*(Ns+1)*(Ns+1);

  unsigned char *tmp_mat = new unsigned char[elemSize];  
  double *tmp_tau = new double[elemSize];
  double *tmp_fib = new double[dof*elemSize];

  // generate filenames & read in the raw arrays first ...
  std::ifstream fin;

  sprintf(filename, "%s.%d.img", problemName, Ns); 
  fin.open(filename, std::ios::binary); fin.read((char *)tmp_mat, elemSize); fin.close();

  // Set Elemental material properties
  PetscScalar ***muArray, ***lambdaArray, ***rhoArray;

  CHKERRQ(DAVecGetArray(da, mu, &muArray));
  CHKERRQ(DAVecGetArray(da, lambda, &lambdaArray));
  CHKERRQ(DAVecGetArray(da, rho, &rhoArray));

  // assign material properties ...
  //       myo,   tissue
  // nu  = 0.49,  0.45
  // E   = 10000, 1000
  // rho = 1.0,   0.1

  // std::cout << "Setting Elemental properties." << std::endl;

  // loop through all elements ...

  nu = 0.45; E = 1000;
  double mmu = E/(2*(1+nu));
  double llam = E*nu/((1+nu)*(1-2*nu));
  nu = 0.45; E = 1000*ctrst;
  double mmu2 = E/(2*(1+nu));
  double llam2 = E*nu/((1+nu)*(1-2*nu));

  if (!rank)
    std::cout << "mu, lam are " << mmu << ", " << llam << std::endl;
  for (int k=z; k<z+zne; k++) {
    for (int j=y; j<y+yne; j++) {
      for (int i=x; i<x+xne; i++) {
        int indx = k*Ns*Ns + j*Ns + i;

        if ( tmp_mat[indx] ) {
          muArray[k][j][i] = mmu2;
          lambdaArray[k][j][i] = llam2;
          rhoArray[k][j][i] = 1.0;
        } else {
          muArray[k][j][i] = mmu;
          lambdaArray[k][j][i] = llam;
          rhoArray[k][j][i] = 1.0;
        }

      } // end i
    } // end j
  } // end k

  // std::cout << "Finished Elemental loop." << std::endl;

  CHKERRQ( DAVecRestoreArray ( da, mu, &muArray ) );
  CHKERRQ( DAVecRestoreArray ( da, lambda, &lambdaArray ) );
  CHKERRQ( DAVecRestoreArray ( da, rho, &rhoArray ) );

  // std::cout << "Finished restoring arrays" << std::endl; 
  // delete temporary buffers

  delete [] tmp_mat;

  // Now set the activation ...
  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));
  // tau = (Vec *) new char*[numSteps+1];

  // std::cout << "Numsteps is " << numSteps << std::endl;
  Vec tauVec, tmpTau;

#ifdef __DEBUG__  
  if (!rank) {
    std::cout << x << ", " << y << ", " << z << " + " << xne << ", " << yne << ", " << zne << std::endl;
  }
#endif  

  PetscScalar ***tauArray;
  // read in the fibers 
  CHKERRQ( DACreateGlobalVector(da3d, &fibers) );
  CHKERRQ( VecSet( fibers, 0.0));


  CHKERRQ(DAVecGetArray(da3d, fibers, &tauArray));

  // std::cout << "Setting force vectors" << std::endl;
  sprintf(filename, "%s.%d.fibers", problemName, Ns);
  // std::cout << "Reading force file " << filename << std::endl;
  std::ifstream fin3(filename, std::ios::binary); fin3.read((char *)tmp_fib, dof*elemSize*sizeof(double)); fin3.close();
  for (int k = z; k < z + zne ; k++) {
    for (int j = y; j < y + yne; j++) {
      for (int i = x; i < x + xne; i++) {
        int indx = dof*((k*(Ns) + j)*(Ns) + i);
        tauArray[k][j][dof*i] = tmp_fib[indx];
        tauArray[k][j][dof*i+1] = tmp_fib[indx+1];
        tauArray[k][j][dof*i+2] = tmp_fib[indx+2];
      }
    }
  }
  CHKERRQ( DAVecRestoreArray ( da3d, fibers, &tauArray ) );

  delete [] tmp_fib;
  // std::cout << "Finished reading fibers" << std::endl;
  // DONE FIBERS

  CHKERRQ( DACreateGlobalVector(da, &tmpTau) );
  double tauNorm;
  for (unsigned int t=0; t<numSteps+1; t++) {
    CHKERRQ( DACreateGlobalVector(da, &tauVec) );
    CHKERRQ( VecSet( tmpTau, 0.0));

    CHKERRQ(DAVecGetArray(da, tmpTau, &tauArray));

    sprintf(filename, "%s.%d.%.3d.fld", problemName, Ns, t);
    std::ifstream fin2(filename); fin2.read((char *)tmp_tau, elemSize*sizeof(double)); fin2.close();

    for (int k = z; k < z + zne ; k++) {
      for (int j = y; j < y + yne; j++) {
        for (int i = x; i < x + xne; i++) {
          int indx = (k*(Ns) + j)*(Ns) + i;
          tauArray[k][j][i] = -10000.0*tmp_tau[indx];
        }
      }
    }
    CHKERRQ( DAVecRestoreArray ( da, tmpTau, &tauArray ) );
    elementToNode(da, tmpTau, tauVec);

    tau.push_back(tauVec);
  }
  CHKERRQ( VecDestroy( tmpTau ) );
  delete [] tmp_tau;

  std::cout << "Finished reading all files" << std::endl;

  // DONE - SET MATERIAL PROPERTIES ...
  // Setup Matrices and Force Vector ...
  Mass->setProblemDimensions(1.0, 1.0, 1.0);
  Mass->setDA(da3d);
  Mass->setDof(dof);
  Mass->setDensity(rho);

  Stiffness->setProblemDimensions(1.0, 1.0, 1.0);
  Stiffness->setDA(da3d);
  Stiffness->setDof(dof);
  Stiffness->setLame(lambda, mu);

  Damping->setAlpha(0.0);
  Damping->setBeta(0.00075);
  Damping->setMassMatrix(Mass);
  Damping->setStiffnessMatrix(Stiffness);
  Damping->setDA(da3d);
  Damping->setDof(dof);

  // Force Vector
  Force->setProblemDimensions(1.0,1.0,1.0);
  Force->setDA(da3d);
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
  ts->storeVec(true);

  ts->setForceVector(Force);

  ts->setInitialDisplacement(initialDisplacement);
  ts->setInitialVelocity(initialVelocity);

  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward
  ts->useMatrixFree(mfree);

  //if (!rank)
  std::cout << RED"Initializing Newmark"NRM << std::endl;
  double itime = MPI_Wtime();

  ts->init(); // initialize IMPORTANT 
  std::cout << "Solving newmark" << std::endl;
  // ts->solve();
  std::cout << "Done solving newmark" << std::endl;


  // Initial guess ...
  Vec alpha, outvec;
  PetscScalar* avec;
  VecCreateSeq(PETSC_COMM_SELF, numParams, &alpha);
  VecDuplicate(alpha, &outvec);
  VecZeroEntries(alpha);
	VecZeroEntries(outvec);

  // Inverse solver set up
  parametricActivationInverse *hyperInv = new parametricActivationInverse;
  // PetscPrintf(0, "Constructed\n");

  hyperInv->setScalarDA(da);

  hyperInv->setBasis(spatialBasis, temporalBasis);
  // PetscPrintf(0, "Set Bases\n");
  hyperInv->setForwardInitialConditions(initialDisplacement, initialVelocity);
  // PetscPrintf(0, "Set Initial\n");
  hyperInv->setTimeStepper(ts);    // set the timestepper
  hyperInv->setInitialGuess(alpha);// set the initial guess 
  hyperInv->setRegularizationParameter(beta); // set the regularization paramter

  // hyperInv->setObservations(solvec); // set the data for the problem 
  hyperInv->init(); // initialize the inverse solver

  hyperInv->solve();

  Vec FinalSolution;
  hyperInv->getCurrentControl(FinalSolution);


  char fname[256];
  sprintf(fname, "%s.soln.%d.%d.raw",problemName, Ns, parFac );

  std::ofstream sol;
  if (!rank) {
    sol.open(fname, std::ios::binary);

    VecGetArray(FinalSolution, &avec);
    sol.write((char *)avec, numParams*sizeof(PetscScalar));
    VecRestoreArray(outvec, &avec);

    sol.close();
  }



  /*
  char fname[256];
  sprintf(fname, "%s.hessFiber.%d.%d.raw",problemName, Ns, parFac );
  std::ofstream hess;

  if (!rank) {
    hess.open(fname, std::ios::binary);
  }

  for (int i=0; i<numParams; i++) {
    std::cout << GRN"Hessian row "RED << i << NRM << std::endl;
    VecZeroEntries(alpha);
    VecGetArray(alpha, &avec);
    avec[i] = 1.0;
    VecRestoreArray(alpha, &avec);


    hyperInv->hessianMatMult(alpha, outvec);
    // VecView(outvec, 0);

    // save the result ...
    if (!rank) {
      VecGetArray(outvec, &avec);
      hess.write((char *)avec, numParams*sizeof(PetscScalar));
      VecRestoreArray(outvec, &avec);
    }
  }
  if (!rank) {
    hess.close();
  }
  */
  
    PetscFinalize();
}


