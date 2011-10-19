static char help[] = "Driver for a linear elastodynamic problem (Hyperbolic)";

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>

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
#include "hyperbolicInverse.h"

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
  double gamma = 0.0; // percent of noise added ...

  // double dtratio = 1.0;
  DA  da;         // Underlying scalar DA - for scalar properties
  DA  da3d;       // Underlying vector DA - for vector properties

  Vec rho;        // density - elemental scalar
  Vec lambda;     // Lame parameter - lambda - elemental scalar
  Vec mu;         // Lame parameter - mu - elemental scalar
  Vec fibers;     // Fiber orientations - nodal vector (3-dof)
  Vec fibersElemental; // for IO. will be destroyed. - elemental vector (3-dof)

  std::vector<Vec> tau;        // the scalar activation - nodal scalar

  // Initial conditions
  Vec initialDisplacement; 
  Vec initialVelocity;

  timeInfo ti;

  // get Ns
  CHKERRQ ( PetscOptionsGetInt(0,"-Ns",&Ns,0) );

  CHKERRQ ( PetscOptionsGetScalar(0,"-t0",&t0,0) );
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

  cardiacDynamic *Force = new cardiacDynamic(feVec::PETSC); // Force Vector

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
  for (int k=z; k<z+zne; k++) {
    for (int j=y; j<y+yne; j++) {
      for (int i=x; i<x+xne; i++) {
        int indx = k*Ns*Ns + j*Ns + i;

        if ( tmp_mat[indx] ) {
          muArray[k][j][i] = 344.82; //3355.7;
          lambdaArray[k][j][i] = 3103.448;// 164429.53;
          rhoArray[k][j][i] = 1.0;
        } else {
          muArray[k][j][i] = 344.82;
          lambdaArray[k][j][i] = 3103.448;
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
          tauArray[k][j][dof*i] = tmp_tau[indx];
          tauArray[k][j][dof*i+1] = tmp_tau[indx+1];
          tauArray[k][j][dof*i+2] = tmp_tau[indx+2];
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
  //if (!rank) {
  //  std::cout << "Finished setting forces" << std::endl;
  // }

  // CHKERRQ( VecDestroy( tmpTau ) );
  delete [] tmp_tau;

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

  ts->setForceVector(Force);

  ts->setInitialDisplacement(initialDisplacement);
  ts->setInitialVelocity(initialVelocity);
  ts->storeVec(true);
  ts->setTimeInfo(&ti);
  ts->setAdjoint(false); // set if adjoint or forward

  ts->init(); // initialize IMPORTANT 
  if (!rank)
    std::cout << RED"Starting Newmark Solve"NRM << std::endl;
  ts->solve();// solve 
  if (!rank)
    std::cout << GRN"Done Newmark"NRM << std::endl;


  std::vector<Vec> solvec = ts->getSolution();

  /* Set very initial guess for the inverse problem*/
  /*
  PetscRandom rctx;
  PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  PetscRandomSetFromOptions(rctx);
  VecSetRandom(guess,rctx);
  VecNorm(guess,NORM_INFINITY,&norm);
  PetscPrintf(0,"guess norm = %g\n",norm);
  */

  double errnorm;
  double exsolnorm;

  Vec guess;
  Vec truth;
  Vec Err;
  concatenateVecs(solvec, guess);
  concatenateVecs(tau, truth);
  
  iC(VecNorm(truth, NORM_2, &exsolnorm));

  /*
  std::cout << "Forward solver solution size is " << solvec.size() << std::endl;
  std::cout << "Forward solver solution norm is " << exsolnorm << std::endl;
*/
  VecZeroEntries(guess);

  // Inverse solver set up
  hyperbolicInverse *hyperInv = new hyperbolicInverse;

  hyperInv->setForwardInitialConditions(initialDisplacement, initialVelocity);
  hyperInv->setTimeStepper(ts);    // set the timestepper
  hyperInv->setInitialGuess(guess);// set the initial guess 
  // hyperInv->setInitialGuess(truth);// set the initial guess 
  hyperInv->setRegularizationParameter(beta); // set the regularization paramter
  hyperInv->setObservations(solvec); // set the data for the problem 
  hyperInv->init(); // initialize the inverse solver
  hyperInv->solve(); // solve
  hyperInv->getCurrentControl(guess); // get the solution 



  // see the error in the solution relative to the actual solution
  VecDuplicate(truth, &Err);
  iC(VecZeroEntries(Err)); 
  iC(VecWAXPY(Err, -1.0, guess, truth));
  iC(VecNorm(Err, NORM_2, &errnorm));
  PetscPrintf(0,"errr in inverse = %g\n", errnorm/exsolnorm);

  PetscFinalize();
}

