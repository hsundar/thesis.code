
static char help[] = "Driver for the parameterized wave equation (Hyperbolic)";

#include "mpi.h"
#include <iostream>
#include <fstream>

#include "petscksp.h"
#include "petscda.h"

#include "timeInfo.h"
#include "feMatrix.h"
#include "feVector.h"
#include "femUtils.h"
#include "timeStepper.h"
#include "newmark.h"
#include "stiffnessMatrix.h"
#include "massMatrix.h"
#include "waveDamping.h"
#include "fdynamicVector.h"
#include "parametricWaveInverse.h"

void getForces(Vec params, std::vector<std::vector<Vec> > forceBasis, std::vector<Vec> &forces);
void getForces(Vec params, std::vector<Vec> &forces, DA da, timeInfo ti, int numParams);

int main(int argc, char **argv)
{    
  PetscInitialize(&argc, &argv, "wave.opt", help);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double startTime, endTime;
  
  int Ns = 32;
  unsigned int dof = 1;

  // double dtratio = 1.0;
  DA  da;         // Underlying DA

  Vec rho;        // density - elemental scalar
  Vec nu;         // Lame parameter - lambda - elemental scalar

  std::vector < std::vector<Vec> > fBasis;        // the scalar activation - nodal scalar
  // std::vector<Vec> truth;        // the ground truth.

  // Initial conditions
  Vec initialDisplacement; 
  Vec initialVelocity;

  timeInfo ti;

  // get Ns
  CHKERRQ ( PetscOptionsGetInt(0,"-Ns", &Ns,0) );

  double t0 = 0.0;
  double dt = 1.0/(Ns);
  double t1 = 1.0;

  double nuVal = 1.0;
  double beta = 0.0001;
  int numParams = 5;

  CHKERRQ ( PetscOptionsGetInt(0,"-nump",&numParams,0) );

  CHKERRQ ( PetscOptionsGetScalar(0,"-t0",&t0,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-t1",&t1,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-dt",&dt,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-nu",&nuVal,0) );
  CHKERRQ ( PetscOptionsGetScalar(0,"-beta",&beta,0) );
  // CHKERRQ ( PetscOptionsGetString(PETSC_NULL, "-pn", problemName, PETSC_MAX_PATH_LEN-1, PETSC_NULL));

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  if (!rank) {
    std::cout << "Problem size is " << Ns+1 << " spatially and NT = " << (int)ceil(1.0/dt) << std::endl << std::endl;
    std::cout << "Number of parameters is " << numParams << std::endl;
  }
  // create DA
  CHKERRQ ( DACreate3d ( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, 
                         Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                         1, 1, 0, 0, 0, &da) );

  massMatrix *Mass = new massMatrix(feMat::PETSC); // Mass Matrix
  stiffnessMatrix *Stiffness = new stiffnessMatrix(feMat::PETSC); // Stiffness matrix
  waveDamping *Damping = new waveDamping(feMat::PETSC); // Damping Matrix

  fdynamicVector *Force = new fdynamicVector(feVec::PETSC); // Force Vector

  // create vectors 
  CHKERRQ( DACreateGlobalVector(da, &rho) );
  CHKERRQ( DACreateGlobalVector(da, &nu) );

  CHKERRQ( DACreateGlobalVector(da, &initialDisplacement) );
  CHKERRQ( DACreateGlobalVector(da, &initialVelocity) );

  // Set initial conditions
  CHKERRQ( VecSet ( initialDisplacement, 0.0) ); 
  CHKERRQ( VecSet ( initialVelocity, 0.0) );

  VecZeroEntries( nu );
  VecZeroEntries( rho );

  CHKERRQ( VecSet ( nu, nuVal) ); 
  CHKERRQ( VecSet ( rho, 1.0) );

  int x, y, z, m, n, p;
  int mx,my,mz;

  CHKERRQ( DAGetCorners(da, &x, &y, &z, &m, &n, &p) ); 
  CHKERRQ( DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ); 

  double acx,acy,acz;
  double hx = 1.0/((double)Ns);

  // allocate for temporary buffers ...
  // unsigned int elemSize = Ns*Ns*Ns;
  // std::cout << "Elem size is " << elemSize << std::endl;
  // unsigned int nodeSize = (Ns+1)*(Ns+1)*(Ns+1);

  // Now set the activation ...
  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));

 // Vec tauVec;

  // PetscScalar ***tauArray;

  unsigned int paramTimeSteps = (unsigned int)(ceil(( (double)(numSteps))/ ((double)(2*numParams)) ));

  /*
  for (int b=0; b<numParams; b++) {
    std::vector<Vec> tau;
    unsigned int tBegin = paramTimeSteps*b;
    unsigned int tEnd   = tBegin + numSteps/2; // paramTimeSteps*(b+2);

    // std::cout << "For param " << b << ": Time step range is " << tBegin << " -> " << tEnd << std::endl; 
    for (unsigned int t=0; t<numSteps+1; t++) {
      double newTime = (dt*(t-tBegin)*numSteps)/((double)(paramTimeSteps));
     // double fff = 0.0;
      CHKERRQ( DACreateGlobalVector(da, &tauVec) );
      CHKERRQ( VecSet( tauVec, 0.0));

      if ( (t>=tBegin) && (t<=tEnd)) {
        CHKERRQ(DAVecGetArray(da, tauVec, &tauArray));
        for (int k = z; k < z + p ; k++) {
          for (int j = y; j < y + n; j++) {
            for (int i = x; i < x + m; i++) {
              acx = (i)*hx; acy = (j)*hx; acz = (k)*hx;
              tauArray[k][j][i] = sin(M_PI*newTime)*cos(2*M_PI*acx)*cos(2*M_PI*acy)*cos(2*M_PI*acz);
            }
          }
        }
        CHKERRQ( DAVecRestoreArray ( da, tauVec, &tauArray ) );
      }
      tau.push_back(tauVec);
    }
    fBasis.push_back(tau);
  }
  */
 // std::cout << "Finished setting basis" << std::endl;

  /*
  // Set initial velocity ...
  CHKERRQ(DAVecGetArray(da, initialVelocity, &solArray));

  for (int k = z; k < z + p ; k++) {
  for (int j = y; j < y + n; j++) {
  for (int i = x; i < x + m; i++) {
  acx = (i)*hx; acy = (j)*hx; acz = (k)*hx;
  solArray[k][j][i] = M_PI*cos(2*M_PI*acx)*cos(2*M_PI*acy)*cos(2*M_PI*acz);
  }
  }
  }
  CHKERRQ( DAVecRestoreArray ( da, initialVelocity, &solArray ) );
  */

  std::vector<Vec> newF;

  Vec alpha;
  PetscScalar *avec;

  VecCreateSeq(PETSC_COMM_SELF, numParams, &alpha);
  /*
  VecCreate(PETSC_COMM_WORLD, &alpha);
  VecSetSizes(alpha, numParams, PETSC_DECIDE);
  VecSetFromOptions(alpha);
  */

  VecGetArray(alpha, &avec);

  for (int j=0; j<numParams; j++)
    avec[j] = 0.5 + 0.5*j;

  VecRestoreArray(alpha, &avec);

  // getForces(alpha, fBasis, newF);
  getForces(alpha, newF, da, ti, numParams);

  // Setup Matrices and Force Vector ...
  Mass->setProblemDimensions(1.0, 1.0, 1.0);
  Mass->setDA(da);
  Mass->setDof(dof);
  Mass->setNuVec(rho);

  Stiffness->setProblemDimensions(1.0, 1.0, 1.0);
  Stiffness->setDA(da);
  Stiffness->setDof(dof);
  Stiffness->setNuVec(nu);

  Damping->setAlpha(0.0);
  Damping->setBeta(0.00075);
  Damping->setMassMatrix(Mass);
  Damping->setStiffnessMatrix(Stiffness);
  Damping->setDA(da);
  Damping->setDof(dof);

  // Force Vector
  Force->setProblemDimensions(1.0,1.0,1.0);
  Force->setDA(da);
  Force->setFDynamic(newF);
  Force->setTimeInfo(&ti);

  // Newmark time stepper ...
  newmark *ts = new newmark; 

  ts->setMassMatrix(Mass);
  ts->setDampingMatrix(Damping);
  ts->setStiffnessMatrix(Stiffness);
  ts->damp(false);
  ts->setTimeFrames(1);
  ts->storeVec(true);
  ts->setAdjoint(false);

  ts->setForceVector(Force);

  ts->setInitialDisplacement(initialDisplacement);
  ts->setInitialVelocity(initialVelocity);

  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward

  ts->init(); // initialize IMPORTANT 
 // if (!rank)
 //   std::cout << RED"Starting initial forward solve"NRM << std::endl;
  ts->solve();// solve 
 // if (!rank)
 // std::cout << GRN"Finished with initial forward solve"NRM << std::endl;

  std::vector<Vec> solvec = ts->getSolution();
  // Now lets check the error ...
  // Vec nr;
  // concatenateVecs(solvec, nr);

  // VecDestroy(nr);
  // VecDestroy(gt);

  // std::cout << std::endl;
  /*************
   *  INVERSE  *
  *************/

  // True solution is tau ... we want to recover it.
  // The observations in this case are, solvec

  /* Set very initial guess for the inverse problem*/

   // Now can clear memory ...

  /*
  for (int i=0; i<newF.size(); i++) {
    if (newF[i] != NULL) {
      VecDestroy(newF[i]);
    }
  }
  newF.clear();

  for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  } 
  solvec.clear();

  ts->destroy();

  VecDestroy(rho);
  VecDestroy(nu);
  VecDestroy(initialDisplacement);
  VecDestroy(initialVelocity);

  VecDestroy(alpha);

  DADestroy(da);

  PetscFinalize();

  return 0;
  */

  Vec gt, nr;
  concatenateVecs(solvec, gt);

  Vec guess;
  VecDuplicate(alpha, &guess);
  VecZeroEntries(guess);
  // VecDuplicate(guess, &Out);
  // VecZeroEntries(Out);

  // double norm;
  /*
     PetscRandom rctx;
     PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
     PetscRandomSetFromOptions(rctx);
     VecSetRandom(guess, rctx);
     VecNorm(guess, NORM_2, &norm);
     PetscPrintf(0, "guess norm = %g\n", norm);
     */


  // double errnorm;
  // double exsolnorm;

  // Inverse solver set up
  // std::cout << RED"Setting up Inverse Solver"NRM << std::endl;
  parametricWaveInverse* hyperInv = new parametricWaveInverse;
  // std::cout << GRN"Finished setting up Inverse Solver"NRM << std::endl;


  hyperInv->setTimeStepper(ts);    // set the timestepper
  hyperInv->setForwardInitialConditions(initialDisplacement, initialVelocity);
  // std::cout << RED"Setting initial guess"NRM << std::endl;
  // hyperInv->setInitialGuess(truth);// set the initial guess 
  hyperInv->setInitialGuess(guess);// set the initial guess 
  // std::cout << GRN"Done setting initial guess"NRM << std::endl;
  hyperInv->setRegularizationParameter(beta); // set the regularization paramter
  hyperInv->setAdjoints(solvec); // set the data for the problem 

  // hyperInv->setForceBasis(fBasis);
  hyperInv->setNumberOfParameter(numParams);

  // std::cout << RED"Initializing Inverse Solver"NRM << std::endl;
  hyperInv->init(); // initialize the inverse solver

 // if (!rank)
 //   std::cout << RED"Starting Inverse Solve"NRM << std::endl;
  startTime = MPI_Wtime();
  hyperInv->solve(); // solve
  endTime = MPI_Wtime();
 // if (!rank)
 //   std::cout << GRN"FINISHED HESSIAN SOLVE"NRM << std::endl;

  hyperInv->getCurrentControl(guess); // get the solution 

  hyperInv->destroy();

  /*
 for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  } 
  solvec.clear();
  */
  // VecView(guess, 0);

  if (!rank)
    std::cout << std::endl << "Error Norms " << std::endl;

  Vec Err;
  double gtNorm, solNorm, errNorm;
  VecDuplicate(guess, &Err);
  VecWAXPY(Err, -1.0, guess, alpha);
  VecNorm(alpha, NORM_2, &gtNorm);
  VecNorm(guess, NORM_2, &solNorm);
  VecNorm(Err, NORM_2, &errNorm);

  if (!rank) {
    std::cout << "The norms are " << gtNorm << ", " << solNorm << ", " << errNorm << std::endl;
    std::cout << "Relative error is " << errNorm/gtNorm << std::endl;
  }
  // Now we shall do another forward solve ...
  getForces(guess, newF, da, ti, numParams);
  Force->setFDynamic(newF);

  ts->setInitialDisplacement(initialDisplacement);
  ts->setInitialVelocity(initialVelocity);
  ts->setAdjoint(false);
  ts->clearMonitor();
  ts->solve();

  std::vector<Vec> solvec2 = ts->getSolution();

  ts->destroy();

  concatenateVecs(solvec2, nr);

   // Now can clear memory ...
  for (int i=0; i<solvec2.size(); i++) {
    if (solvec2[i] != NULL) {
      VecDestroy(solvec2[i]);
    }
  }
  solvec2.clear();

   // Now can clear memory ...
  for (int i=0; i<newF.size(); i++) {
    if (newF[i] != NULL) {
      VecDestroy(newF[i]);
    }
  }
  newF.clear();
/*
  for (unsigned int i=0; i<truth.size(); i++) {
    VecNorm(truth[i], NORM_2, &gtNorm);
    VecNorm(solvec[i], NORM_2, &solNorm);
    VecAXPY(solvec[i], -1.0, truth[i]);
    VecNorm(solvec[i], NORM_2, &errNorm);
    PetscPrintf(0, "Ground truth at timestep %d is %g, %g, %g\n", i, gtNorm, solNorm, errNorm);
    // PetscPrintf(0, "Relative Error at timestep %d is %g\n", i, errNorm/gtNorm);
  }
  */
  VecNorm(gt, NORM_2, &gtNorm);
  VecAXPY(nr, -1.0, gt);
  VecNorm(nr, NORM_2, &errNorm);

  if (!rank)
    std::cout <<  "Total Relative error on state is " << errNorm/gtNorm << std::endl;
  
  if (!rank)
    std::cout << "Wall time is " << endTime - startTime << std::endl;


  VecDestroy(gt);
  VecDestroy(nr);
  VecDestroy(Err);
  VecDestroy(alpha);
  VecDestroy(guess);

  VecDestroy(rho);
  VecDestroy(nu);
  VecDestroy(initialDisplacement);
  VecDestroy(initialVelocity);

  DADestroy(da);

  PetscFinalize();
}

void getForces(Vec params, std::vector<Vec> &forces, DA da, timeInfo ti, int numParams) {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif
  // Generate the force vector based on the current parameters ...
  // F = Sum { B_i a_i}
  
  PetscScalar * pVec;
  VecGetArray(params, &pVec);
  // Clear the Forces
  for (unsigned int i=0; i<forces.size(); i++) {
    if (forces[i] != NULL) {
      VecDestroy(forces[i]);
    }
  }
  forces.clear();

  unsigned int numSteps = (unsigned int)(ceil(( ti.stop - ti.start)/ti.step));
  // create and initialize to 0
  for (unsigned int i=0; i<numSteps+1; i++) {
    Vec tmp;
    DACreateGlobalVector(da, &tmp);
    VecZeroEntries(tmp);
    forces.push_back(tmp);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  PetscScalar ***tauArray;

  unsigned int paramTimeSteps = (unsigned int)(ceil(( (double)(numSteps))/ ((double)(2*numParams)) ));
  double acx,acy,acz;

  int x, y, z, m, n, p;
  int mx,my,mz;

  DAGetCorners(da, &x, &y, &z, &m, &n, &p);
  DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0);

  double hx = 1.0/(mx-1.0);

  for (int b=0; b<numParams; b++) {
    std::vector<Vec> tau;
    unsigned int tBegin = paramTimeSteps*b;
    unsigned int tEnd   = tBegin + numSteps/2; // paramTimeSteps*(b+2);

    // std::cout << "For param " << b << ": Time step range is " << tBegin << " -> " << tEnd << std::endl; 
    for (unsigned int t=0; t<numSteps+1; t++) {
      double newTime = (ti.step*(t-tBegin)*numSteps)/((double)(paramTimeSteps));
      
      if ( (t>=tBegin) && (t<=tEnd)) {
        DAVecGetArray(da, forces[t], &tauArray);
        for (int k = z; k < z + p ; k++) {
          for (int j = y; j < y + n; j++) {
            for (int i = x; i < x + m; i++) {
              acx = (i)*hx; acy = (j)*hx; acz = (k)*hx;
              tauArray[k][j][i] += pVec[b]*sin(M_PI*newTime)*cos(2*M_PI*acx)*cos(2*M_PI*acy)*cos(2*M_PI*acz);
            }
          }
        }
        DAVecRestoreArray ( da, forces[t], &tauArray ) ;
      }
    }
  }
  VecRestoreArray(params, &pVec);

#ifdef __DEBUG__
  // Get the norms of the forces ... just to be safe ..
  double fNorm1, fNorm2;

  for (unsigned int i=0; i<forces.size(); i++) {
    VecNorm(forces[i], NORM_INFINITY, &fNorm1);
    VecNorm(forces[i], NORM_2, &fNorm2);
    PetscPrintf(0, "Force Norms at timestep %d are %g and %g\n", i, fNorm1, fNorm2);
  }
#endif
  
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}

// When the discrete forces are given ...
void getForces(Vec params, std::vector<std::vector<Vec> > forceBasis, std::vector<Vec> &forces) {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif
  // Generate the force vector based on the current parameters ...
  // F = Sum { B_i a_i}
  
  PetscScalar * pVec;
  VecGetArray(params, &pVec);
  // Clear the Forces
  for (unsigned int i=0; i<forces.size(); i++) {
    if (forces[i] != NULL) {
      VecDestroy(forces[i]);
    }
  }
  forces.clear();

  // create and initialize to alpha_0 * Basis_0
  for (unsigned int i=0; i<forceBasis[0].size(); i++) {
    Vec tmp;
    VecDuplicate(forceBasis[0][i], &tmp);
    VecZeroEntries(tmp);
    VecAXPY(tmp, pVec[0], forceBasis[0][i]);
    forces.push_back(tmp);
  }

  for (unsigned int i=1; i<forceBasis.size(); i++) {
    for (unsigned int j=0; j<forceBasis[i].size(); j++) {
      VecAXPY(forces[j], pVec[i], forceBasis[i][j]);
    }
  }
  VecRestoreArray(params, &pVec);

#ifdef __DEBUG__
  // Get the norms of the forces ... just to be safe ..
  double fNorm1, fNorm2;

  for (unsigned int i=0; i<forces.size(); i++) {
    VecNorm(forces[i], NORM_INFINITY, &fNorm1);
    VecNorm(forces[i], NORM_2, &fNorm2);
    PetscPrintf(0, "Force Norms at timestep %d are %g and %g\n", i, fNorm1, fNorm2);
  }
#endif
  
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}


