static char help[] = "Generate the parametrized linear elastodynamic problem (Hyperbolic)";

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
#include "cardiacDynamic.h"
#include "parametricElasInverse.h"
#include "radialBasis.h"
#include "bSplineBasis.h"

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

void getForces(Vec params, std::vector<Vec> &forces, DA da, timeInfo *ti, std::vector<radialBasis> gBasis, bSplineBasis bBasis);

int main(int argc, char **argv) {
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
  DA  da;					// Underlying scalar DA - for scalar properties
  DA  da3d;				// Underlying vector DA - for vector properties

  Vec rho;				// density - elemental scalar
  Vec lambda;			// Lame parameter - lambda - elemental scalar
  Vec mu;					// Lame parameter - mu - elemental scalar
  Vec fibers;			// Fiber orientations - nodal vector (3-dof)
  Vec fibersElemental; // for IO. will be destroyed. - elemental vector (3-dof)

  std::vector<Vec> tau;				 // the scalar activation - nodal scalar

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
    mfree	= false;

  double ctrst = 1.0;
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

  // Time info for timestepping
  ti.start = t0;
  ti.stop  = t1;
  ti.step  = dt;

  if (!rank) {
    std::cout << "Grid size is " << Ns+1 << " and NT is " << (int)ceil(1.0/dt) << std::endl;
  }

  // create DA
  CHKERRQ ( DACreate3d ( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, 
        Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
        1, 1, 0, 0, 0, &da) );
  CHKERRQ ( DACreate3d ( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, 
        Ns+1, Ns+1, Ns+1, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
        dof, 1, 0, 0, 0, &da3d) );

  elasMass *Mass = new elasMass(feMat::PETSC); // Mass Matrix
  elasStiffness *Stiffness = new elasStiffness(feMat::PETSC);	// Stiffness matrix
  raleighDamping *Damping = new raleighDamping(feMat::PETSC);	// Damping Matrix

  cardiacDynamic *Force = new cardiacDynamic(feVec::PETSC);	// Force Vector

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

  int parFac = 2;
  int numParams = 120;

  CHKERRQ ( PetscOptionsGetInt(0,"-pFac", &parFac,0) );

  numParams = parFac*parFac*parFac*5*3;
  if (!rank)
    std::cout << "Total number of unknowns is " << numParams << std::endl;

  // Generate the basis ...
  std::vector < radialBasis > spatialBasis;
  bSplineBasis temporalBasis(3, 5);	// this creates and sets up the basis ...
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
  // std::cout << "Elem size is " << elemSize << std::endl;
  unsigned int nodeSize = (Ns+1)*(Ns+1)*(Ns+1);

  unsigned char *tmp_mat = new unsigned char[elemSize];  
  double *tmp_tau = new double[dof*elemSize];

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
  nu = 0.45; E = 1000;
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

      }	// end i
    }	// end j
  }	// end k

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

  CHKERRQ( DACreateGlobalVector(da3d, &tmpTau) );

#ifdef __DEBUG__  
  if (!rank) {
    std::cout << x << ", " << y << ", " << z << " + " << xne << ", " << yne << ", " << zne << std::endl;
  }
#endif  

  PetscScalar ***tauArray;
  double tauNorm;
  for (unsigned int t=0; t<numSteps+1; t++) {
    CHKERRQ( DACreateGlobalVector(da3d, &tauVec) );
    CHKERRQ( VecSet( tmpTau, 0.0));

    CHKERRQ(DAVecGetArray(da3d, tmpTau, &tauArray));

    // std::cout << "Setting force vectors" << std::endl;
    sprintf(filename, "%s.%d.%.3d.fld", problemName, Ns, t);
    // std::cout << "Reading force file " << filename << std::endl;
    fin.open(filename); fin.read((char *)tmp_tau, dof*elemSize*sizeof(double)); fin.close();
    for (int k = z; k < z + zne ; k++) {
      for (int j = y; j < y + yne; j++) {
        for (int i = x; i < x + xne; i++) {
          int indx = dof*(k*Ns*Ns + j*Ns + i);
          tauArray[k][j][dof*i] = 30000*tmp_tau[indx];
          tauArray[k][j][dof*i+1] = 30000*tmp_tau[indx+1];
          tauArray[k][j][dof*i+2] = 30000*tmp_tau[indx+2];
        }
      }
    }
    // std::cout << CYN"\tFinished elemental loop"NRM << std::endl;
    CHKERRQ( DAVecRestoreArray ( da3d, tmpTau, &tauArray ) );
    // std::cout << "Converting to Nodal Vector" << std::endl;
    // VecNorm(tmpTau, NORM_2, &tauNorm);
    // tauNorm = tauNorm/pow(Ns,1.5);
    // std::cout << "Elemental Norm is " << tauNorm << std::endl;
    // std::cout << rank << " Converting to Nodal" << std::endl;
    elementToNode(da3d, tmpTau, tauVec);
    /*
       VecNorm(tauVec, NORM_2, &tauNorm);
       tauNorm = tauNorm/pow(Ns,1.5);
       std::cout << "Nodal Norm is " << tauNorm << std::endl;
       */

    // std::cout << rank << " Done converting to Nodal Vector" << std::endl;
    tau.push_back(tauVec);
  }
  // if (!rank) {
  //   std::cout << "Finished setting forces" << std::endl;
  // }

  // CHKERRQ( VecDestroy( tmpTau ) );
  delete [] tmp_tau;

  // DONE - SET MATERIAL PROPERTIES ...
  // std::cout << "Setting material properties" << std::endl;
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
  // Force->setActivationVec(tau);
  // Force->setFiberOrientations(fibers);
  Force->setFDynamic(tau);
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

  ts->init();	// initialize IMPORTANT 

  // PetscFinalize();
  Vec alpha, outvec;
  PetscScalar *avec;

  VecCreateSeq(PETSC_COMM_SELF, numParams, &alpha);

  VecGetArray(alpha, &avec);

  for (int j=0; j<numParams; j++)
    avec[j] = 0.0;

  avec[1] = 1.0;

  VecRestoreArray(alpha, &avec);

  VecDuplicate(alpha, &outvec);

  // Inverse solver set up
  parametricElasInverse *hyperInv = new parametricElasInverse;

  hyperInv->setBasis(spatialBasis, temporalBasis);

  hyperInv->setForwardInitialConditions(initialDisplacement, initialVelocity);

  hyperInv->setTimeStepper(ts);		 // set the timestepper

  hyperInv->setInitialGuess(alpha);// set the initial guess 

  hyperInv->setRegularizationParameter(beta);	// set the regularization paramter

  // hyperInv->setObservations(solvec); // set the data for the problem 

  hyperInv->init();	// initialize the inverse solver

  // Compute Hessian and write to file ...
  // hyperInv->computeHessian("test.m");

  // std::ofstream hess("hessian32.raw", std::ios::binary);
  char fname[256];
  sprintf(fname, "%s.hessFull.%d.%d.raw",problemName, Ns, parFac );

  std::ofstream hess;
  if (!rank)
    hess.open(fname, std::ios::binary);

  for (int i=0; i<numParams; i++) {
    if (!rank)
      std::cout << GRN"Hessian row "RED << i << NRM << std::endl;
    VecZeroEntries(alpha);
    VecGetArray(alpha, &avec);
    // for (int j=0; j<numParams; j++)
    //	avec[j] = 0.0;
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
  if (!rank)
    hess.close();

  PetscFinalize();
}

