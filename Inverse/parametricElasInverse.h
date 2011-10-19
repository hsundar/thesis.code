/**
 * @file parametricElasInverse.h
 * @brief Main class for an inverse scalar Hyperbolic problem
 * @author Hari Sundar 
 * @date   12/17/07 
 * 
 * Main class for inverse scalar Hyperbolic problem where the 
 * force is parametrised. 
 * 
 **/

#ifndef _PARAMETRIC_ELAS_INVERSE_H_
#define _PARAMETRIC_ELAS_INVERSE_H_

#include <vector> 
#include "inverseSolver.h"
#include "femUtils.h"
#include "radialBasis.h"
#include "bSplineBasis.h"

class parametricElasInverse : public inverseSolver {

public:

  parametricElasInverse() {
  }

  ~parametricElasInverse() {
  }

  virtual int destroy() {
    // Allocate memory for working vectors
    CHKERRQ(VecDestroy(m_vecCurrentControl));
    CHKERRQ(VecDestroy(m_vecControlStep));
    CHKERRQ(VecDestroy(m_vecReducedGradient));

    CHKERRQ(MatDestroy(m_matReducedHessian));
    CHKERRQ(KSPDestroy(m_ksp));

    return true;
  }

  virtual int init();

  virtual int solve();

  virtual void hessianMatMult(Vec In, Vec Out);

  virtual bool setReducedGradient();

  void setForwardInitialConditions(Vec initDisp, Vec initVel) {
    VecDuplicate(initDisp, &m_vecForwardInitialDisplacement);
    VecDuplicate(initVel, &m_vecForwardInitialVelocity);
    VecCopy(initDisp, m_vecForwardInitialDisplacement);
    VecCopy(initVel, m_vecForwardInitialVelocity);
  }

  // MG not yet implemented
  virtual void mghessianMatMult(DA _da, Vec _in, Vec _out) {
  }

  PetscErrorCode setObservations(std::vector<Vec> obs) {
    m_vecObservations = obs;
    return(0);
  }

  void setBasis(std::vector<radialBasis> rb, bSplineBasis bsb) {
    m_radialBasis = rb;
    m_bsplineBasis = bsb;
  }

  // Functions to hanlde the full / parametrized representations
  void getForces(Vec params, std::vector<Vec> &forces);
  void getParams(std::vector<Vec> forces, Vec params);

  void computeHessian(char* fname) {
    std::cout << "Entering computeHessian" << std::endl;
    Mat Hessian;
    KSPComputeExplicitOperator(m_ksp, &Hessian);

    PetscViewer viewer;
    char viewername[100];
    sprintf(viewername, fname);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,viewername, &viewer);
    PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    MatView(Hessian, viewer);
    PetscViewerDestroy(viewer);
    std::cout << "Finished writing Hessian" << std::endl;

  }

protected:
  // In case of time dependent problems, it's faster to save the timesteps separately.
  std::vector<Vec> m_vecObservations;

  // The parametrized force vectors ...
  std::vector<radialBasis> m_radialBasis;
  bSplineBasis             m_bsplineBasis;
  int m_iNumKnots;

  // Inital Vel and Disp for fwd. problem
  Vec m_vecForwardInitialDisplacement;
  Vec m_vecForwardInitialVelocity;



};

/**
 *	@brief The initialization function where the reduced Hessian matrix, KSP context are created
 * @return 0 if successful, 1 otherwise
 * Matrix of shell type is created which does only a matvec
 **/


#undef __FUNCT__
#define __FUNCT__ "pElasInv_Init"
int parametricElasInverse::init() {
#ifdef __DEBUG__  
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif

  int matsize;
  int ierr;

  // Get reduced Hessian matrix size
  ierr = VecDuplicate(m_vecInitialControl, &m_vecCurrentControl); CHKERRQ(ierr);
  ierr = VecDuplicate(m_vecInitialControl, &m_vecReducedGradient); CHKERRQ(ierr);
  ierr = VecDuplicate(m_vecInitialControl, &m_vecControlStep); CHKERRQ(ierr);
  ierr = VecGetLocalSize(m_vecInitialControl, &matsize); CHKERRQ(ierr);

  // std::cout << "Matsize is " << matsize << std::endl;
  // Create Reduced Hessian
  ierr = MatCreateShell(PETSC_COMM_SELF, matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE, this, &m_matReducedHessian); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_matReducedHessian, MATOP_MULT,(void(*)(void))(MatMult)); CHKERRQ(ierr);

  // Create a KSP context to solve for the reduced Hessian system
  ierr = KSPCreate(PETSC_COMM_SELF, &m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(m_ksp,"inv_"); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_ksp, m_matReducedHessian, m_matReducedHessian, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetType(m_ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
  return(0);
}


#undef __FUNCT__
#define __FUNCT__ "pElasInv_Solve"
int parametricElasInverse::solve() {
  int ierr;

  // Set the initial guess to the current control
  ierr = VecCopy(m_vecInitialControl, m_vecCurrentControl); CHKERRQ(ierr);

#ifdef __DEBUG__
  VecNorm(m_vecCurrentControl,NORM_2,&norm);
  PetscPrintf(0,"norm of initial guess = %g\n",norm);
#endif

  // initiate the step to zero
  ierr = VecZeroEntries(m_vecControlStep); CHKERRQ(ierr);

  setReducedGradient();

  // Solve for the step using the reduced Hessian
  ierr = KSPSolve(m_ksp, m_vecReducedGradient, m_vecControlStep); CHKERRQ(ierr);

  PetscReal rnorm;
  PetscInt its;
  KSPGetResidualNorm(m_ksp, &rnorm);
  KSPGetIterationNumber(m_ksp,&its);
  // Print the final residual ...
  PetscPrintf(0, "Final residual norm is %g\n", rnorm);
  PetscPrintf(0, "Total number of iterations: %d\n", its);

#ifdef __DEBUG__
  Mat Hessian; 
  ierr = KSPComputeExplicitOperator(m_ksp,&Hessian); CHKERRQ(ierr);
  std::cout << RED"HESSIAN"NRM << std::endl;
  MatView(Hessian, 0);
#endif

  // new control = old control - step
  ierr = VecAXPY(m_vecCurrentControl, -1.0, m_vecControlStep); CHKERRQ(ierr);

#ifdef __DEBUG__
  VecNorm(m_vecCurrentControl,NORM_INFINITY,&norm);
  PetscPrintf(0,"norm of solution = %g\n",norm);
#endif
  return(0);
}


/**
 *	@brief set the reduced Gradient for the current problem, this function does a forward solve and then -C^T J^-T( y - y*)
 * 
 * Two full solves (one forward + one adjoint)
 **/
#undef __FUNCT__
#define __FUNCT__ "pElasInv_setReducedGradient"
bool parametricElasInverse::setReducedGradient() {

  // std::cout << "entering set RG" << std::endl;
  // Initiate the reduced Gradient to zero
  VecZeroEntries(m_vecReducedGradient);

  // Fdynamic which is the control 
  cardiacDynamic *Fdynamic = new cardiacDynamic(feVec::PETSC);

  // Solver 
  newmark *ts = new newmark;



  // Forward solve
  ts = (newmark *)m_ts;
  ts->setAdjoint(false);

  ts->setInitialDisplacement(m_vecForwardInitialDisplacement);
  ts->setInitialVelocity(m_vecForwardInitialVelocity);

  ts->setTimeFrames(1);
  ts->storeVec(true);

  // Get the force from the timestepper
  Fdynamic = (cardiacDynamic *)ts->getForce();

  // Set the force to the current control
  std::vector<Vec> currControl;
  // splitVec(m_vecCurrentControl, currControl, NT+1);
  // std::cout << "Getting forces" << std::endl;
  getForces(m_vecCurrentControl, currControl);
  Fdynamic->setFDynamic(currControl);

  // set the force in the timestepper
  ts->setForceVector(Fdynamic);

  // clear the monitor
  ts->clearMonitor();
  // ts->init();

  // solve the forward problem to get the state variable @ the current control
  ts->solve();
  // get the solution
  std::vector<Vec> solvec;
  solvec = ts->getSolution();

  // std::cout << "Finished forward solve" << std::endl;
  // Now can clear memory ...
  for (unsigned int i=0; i<currControl.size(); i++) {
    if (currControl[i] != NULL) {
      VecDestroy(currControl[i]);
    }
  }
  currControl.clear();

#ifdef __DEBUG__
  PetscPrintf(0,"size of solvec is %d\n", solvec.size());
  VecNorm(solvec[0],NORM_INFINITY,&norm);
  PetscPrintf(0,"norm of state in reduced gradient = %f\n", norm);
#endif

  // Set the right hand side of the adjoint, (state - data)
  for (unsigned int i=0; i<solvec.size(); i++) {
    VecAYPX(solvec[i], -1.0, m_vecObservations[i]);  
    if (m_bUsePartialObservations) {
      VecPointwiseMult(solvec[i], solvec[i], m_vecPartialObservations);
    }
  }

  // Now can clear memory ...
  for (unsigned int i=0; i<m_vecObservations.size(); i++) {
    if (m_vecObservations[i] != NULL) {
      VecDestroy(m_vecObservations[i]);
    }
  }
  m_vecObservations.clear();

  // set the Fstatic again for the adjoint right hand side
  Fdynamic->setFDynamic(solvec);

  // set the Fstatic
  ts->setForceVector(Fdynamic);

  // set the adjoint flag.. here it does not matter 
  ts->setAdjoint(true);

  // set initial conditions for adjoint ...
  Vec initD, initV;
  VecDuplicate(m_vecForwardInitialDisplacement, &initD);
  VecDuplicate(m_vecForwardInitialDisplacement, &initV);
  VecZeroEntries(initD);
  VecZeroEntries(initV);
  ts->setInitialDisplacement(initD); 
  ts->setInitialVelocity(initV); 

  // Clear monitor
  ts->clearMonitor();
  // ts->init();
  // solve the adjoint problem.. since this is a
  ts->solve();

  // clear memory ...
  for (unsigned int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }
  solvec.clear();

  // get the solution
  solvec = ts->getSolution();

  PetscPrintf(0, "Getting params from forces\n");
  getParams(solvec, m_vecReducedGradient);
  PetscPrintf(0, "Got params from forces\n");

  // clear memory ...
  for (unsigned int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }
  solvec.clear();

  // scale the reduced Gradient
  VecScale(m_vecReducedGradient,-1.0);

  // add the contribution of the regularization parameter
  VecAXPY(m_vecReducedGradient, m_beta, m_vecCurrentControl);

#ifdef __DEBUG__
  VecNorm(m_vecReducedGradient,NORM_2,&norm);
  PetscPrintf(0,"norm of reduced Gradient = %f\n", norm);
#endif
  VecDestroy(initD);
  VecDestroy(initV);

  PetscPrintf(0, "Finished setting reduced gradient\n");
  return true;
}

/**
 *	@brief set the reduced Hessian matrix vector product, this function does one forward solve and one adjoint solve
 *
 * Two full solves : one forward + one adjoint
 **/
#undef __FUNCT__
#define __FUNCT__ "pElasInv_hessianMatMult"
void parametricElasInverse::hessianMatMult(Vec In, Vec Out) {
  // std::cout << "Entering HMM" << std::endl;
  VecZeroEntries(Out);

  newmark *ts = new newmark;

  // Fdynamic which is the control 
  cardiacDynamic *Fdynamic = new cardiacDynamic(feVec::PETSC);

  // Forward solve
  ts = (newmark *)m_ts;

  ts->setInitialDisplacement(m_vecForwardInitialDisplacement);
  ts->setInitialVelocity(m_vecForwardInitialVelocity);

  ts->setTimeFrames(1);
  ts->storeVec(true);
  ts->setAdjoint(false);

  // Get the force from the timestepper
  Fdynamic = (cardiacDynamic *)ts->getForce();

  // Set the force to In
  std::vector<Vec> currControl;
  // splitVec(In, currControl, NT+1);
  getForces(In, currControl);
  Fdynamic->setFDynamic(currControl);

  ts->setForceVector(Fdynamic);
  ts->clearMonitor();
  // ts->init();
  // std::cout << "Starting hessian-forward solve" << std::endl;
  ts->solve();

  // state variable
  std::vector<Vec> solvec;
  solvec = ts->getSolution();

  // Now can clear memory ...
  for (unsigned int i=0; i<currControl.size(); i++) {
    if (currControl[i] != NULL) {
      VecDestroy(currControl[i]);
    }
  }
  currControl.clear();

  // Scale to the set the right hand side of adjoint
  for (unsigned int i=0; i<solvec.size(); i++) {
    VecScale(solvec[i],-1.0);
    if (m_bUsePartialObservations) {
      VecPointwiseMult(solvec[i], solvec[i], m_vecPartialObservations);
    }
  }

  // Adjoint solve steps
  Fdynamic->setFDynamic(solvec);

  // set the adjoint right hand side
  ts->setForceVector(Fdynamic);
  ts->clearMonitor();

  // set the adjoint flag
  ts->setAdjoint(true);

  // set initial conditions for adjoint ...
  Vec initD, initV;
  VecDuplicate(m_vecForwardInitialDisplacement, &initD);
  VecDuplicate(m_vecForwardInitialDisplacement, &initV);
  VecZeroEntries(initD);
  VecZeroEntries(initV);
  ts->setInitialDisplacement(initD); 
  ts->setInitialVelocity(initV); 

  // adjoint solve
  // std::cout << "Starting hessian-adjoint solve" << std::endl;
  ts->solve();

  // clear memory ...
  for (unsigned int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }
  solvec.clear();

  // Adjoint variable
  solvec = ts->getSolution();

  // Now concatenate the solution ...
  getParams(solvec, Out);

  // clear memory ...
  for (unsigned int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }
  solvec.clear();

  // scaling the adjoint variable
  VecScale(Out, -1.0);

  // add the contribution of the regularization parameter
  VecAXPY(Out, m_beta, In);
  // std::cout << GRN"Leaving "NRM << __func__ << std::endl;
  VecDestroy(initD);
  VecDestroy(initV);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Functions to hanlde the full / parametrized representations
#undef __FUNCT__
#define __FUNCT__ "pElasInv_getForces"

void parametricElasInverse::getForces(Vec params, std::vector<Vec> &forces) {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif
  // Generate the activation vector based on the current parameters ...

  PetscScalar * pVec;
  VecGetArray(params, &pVec);

  // Clear the Forces
  for (unsigned int i=0; i<forces.size(); i++) {
    if (forces[i] != NULL) {
      VecDestroy(forces[i]);
    }
  }
  forces.clear();

  double gx;

  // GET DA type ...
  int daType = m_ts->getMass()->getDAtype();

  if ( !daType ) { // type = PETSC
    DA da = m_ts->getMass()->getDA();
    timeInfo *ti = m_ts->getTimeInfo();

    unsigned int numSteps = (unsigned int)(ceil(( ti->stop - ti->start)/ti->step));
    PetscScalar ***tauArray; 

    // create and initialize to 0
    for (unsigned int i=0; i<numSteps+1; i++) {
      Vec tmp;
      DACreateGlobalVector(da, &tmp);
      VecZeroEntries(tmp);
      forces.push_back(tmp);
    }

    int x, y, z, m, n, p;
    int mx,my,mz;

    DAGetCorners(da, &x, &y, &z, &m, &n, &p);
    DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0);

    double hx = 1.0/(mx-1.0);

    unsigned int knotsize = m_bsplineBasis.getNumKnots();
    unsigned int bsize = 3*knotsize;
    double *currBasis = new double[knotsize];

    // std::cout << numSteps << " " << x << " " << y << " " << z << " " << m << " " << n << " " << p << std::endl;
    double currTime = ti->start;
    currTime = ti->start;
    for (unsigned int t=0; t<numSteps+1; t++) {
      m_bsplineBasis.basis(currTime, currBasis);
      DAVecGetArray(da, forces[t], &tauArray);
      for (int k = z; k < z + p ; k++) {
        for (int j = y; j < y + n; j++) {
          for (int i = x; i < x + m; i++) {
            gx = 0.0;
            Point px(i,j,k);
            px *= hx;
            // determine the factor from spatial basis
            for (unsigned int b=0; b<m_radialBasis.size(); b++) {
              double val = m_radialBasis[b].getValue(px);
              gx = val;

              // time loop 
              for (unsigned int r=0; r<knotsize; r++) {
                tauArray[k][j][3*i]    += gx*currBasis[r]*pVec[b*bsize + 3*r ];
                tauArray[k][j][3*i+1]  += gx*currBasis[r]*pVec[b*bsize + 3*r + 1];
                tauArray[k][j][3*i+2]  += gx*currBasis[r]*pVec[b*bsize + 3*r + 2];
              }
            }
          }
        }
      }
      DAVecRestoreArray ( da, forces[t], &tauArray ) ;


      currTime += ti->step;
    }

    delete [] currBasis;

  } else { // otk ...
    ot::DA *da = m_ts->getMass()->getOctDA();

    timeInfo *ti = m_ts->getTimeInfo();

    unsigned int numSteps = (unsigned int)(ceil(( ti->stop - ti->start)/ti->step));
    PetscScalar *tauArray; 

    // create and initialize to 0
    for (unsigned int i=0; i<numSteps+1; i++) {
      Vec tmp;
      da->createVector(tmp, false, true, 3);
      VecZeroEntries(tmp);
      forces.push_back(tmp);
    }

    unsigned int knotsize = m_bsplineBasis.getNumKnots();
    unsigned int bsize = 3*knotsize;
    double *currBasis = new double[knotsize];

    double currTime = ti->start;
    currTime = ti->start;
    for (unsigned int t=0; t<numSteps+1; t++) {
      m_bsplineBasis.basis(currTime, currBasis);

      // DAVecGetArray(da, forces[t], &tauArray);
      da->vecGetBuffer(forces[t], tauArray, false, true, false, 3);

      unsigned int maxD = da->getMaxDepth();
      unsigned int balOctmaxD = maxD - 1;

      for ( da->init<ot::DA::ALL>(), da->init<ot::DA::WRITABLE>(); da->curr() < da->end<ot::DA::ALL>(); da->next<ot::DA::ALL>()) {

        unsigned int i = da->curr();
        Point pt;
        pt = da->getCurrentOffset();
        unsigned levelhere = da->getLevel(da->curr()) - 1;
        double hxOct = (double)((double)(1<<(balOctmaxD - levelhere))/(double)(1 << balOctmaxD));

        double x = (double)(pt.xint())/((double)(1<<(maxD-1)));
        double y = (double)(pt.yint())/((double)(1<<(maxD-1)));
        double z = (double)(pt.zint())/((double)(1<<(maxD-1)));

        Point px(x,y,z);

        for (unsigned int b=0; b<m_radialBasis.size(); b++) {
          double val = m_radialBasis[b].getValue(px);
          gx = val;

          // time loop 
          for (unsigned int r=0; r<knotsize; r++) {
            tauArray[3*i]    += gx*currBasis[r]*pVec[b*bsize + 3*r ];
            tauArray[3*i+1]  += gx*currBasis[r]*pVec[b*bsize + 3*r + 1];
            tauArray[3*i+2]  += gx*currBasis[r]*pVec[b*bsize + 3*r + 2];
          }
        }
      }


      da->vecRestoreBuffer(forces[t], tauArray, false, true, false, 3);

      currTime += ti->step;
    }

    delete [] currBasis;

  }// end else OTK

  VecRestoreArray(params, &pVec);

#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}




#undef __FUNCT__
#define __FUNCT__ "pElasInv_getParams"

void parametricElasInverse::getParams(std::vector<Vec> forces, Vec params) {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif

  int daType = m_ts->getMass()->getDAtype();
  timeInfo *ti = m_ts->getTimeInfo();

  unsigned int numSteps = (unsigned int)(ceil(( ti->stop - ti->start)/ti->step));

  int sz;

  VecZeroEntries(params);
  PetscScalar * pVec;
  VecGetLocalSize(params, &sz);
  PetscScalar *sendVec = new PetscScalar[sz];
  VecGetArray(params, &pVec);

  for (int i=0; i<sz; i++) {
    sendVec[i] = 0.0;
  }

  // number of bSpline bases is dof*knotsize. 
  int knotsize = m_bsplineBasis.getNumKnots();
  int bsize = 3*knotsize;
  double *currBasis = new double[knotsize];

  std::vector<Vec> ibldt;
  for (int i=0; i<knotsize; i++) {
    Vec tmp;
    VecDuplicate(forces[i], &tmp);
    VecZeroEntries(tmp);
    ibldt.push_back(tmp);
  }

  double gx;
  double dt = ti->step;

  // first the time integration of bSpline Basis with lambda
  double currTime = ti->start;
  double lnorm;
  for (unsigned int t=0; t<numSteps+1; t++) {
    m_bsplineBasis.basis(currTime, currBasis);
    for (unsigned int b=0; b<knotsize; b++) {
      double ifac = currBasis[b]*dt; 
      VecAXPY(ibldt[b], ifac, forces[t]);
    }
    currTime += ti->step;
  }


  if ( !daType ) { // PetSc
    PetscScalar ***bl; 

    DA da = m_ts->getMass()->getDA();

    int x, y, z, m, n, p;
    int mx,my,mz;

    DAGetCorners(da, &x, &y, &z, &m, &n, &p);
    DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0);

    double hx = 1.0/(mx-1.0);

    // Now the spatial reduction using the gaussian basis.
    for (unsigned int g=0; g<m_radialBasis.size(); g++) {
      for (int b=0; b<knotsize; b++) {
        DAVecGetArray(da, ibldt[b], &bl);
        for (int k=z; k<z+p; k++) {
          for (int j=y; j<y+n; j++) {
            for (int i=x; i<x+m; i++) {
              Point px(i,j,k);
              px *= hx;
              gx = m_radialBasis[g].getValue(px);

              sendVec[g*bsize + 3*b]     += bl[k][j][3*i]*gx;
              sendVec[g*bsize + 3*b + 1] += bl[k][j][3*i+1]*gx;
              sendVec[g*bsize + 3*b + 2] += bl[k][j][3*i+2]*gx;
            } // i
          } // j
        } // k
        DAVecRestoreArray ( da, ibldt[b], &bl ) ;
      } // b
    } // g
  } else { // OTK
    PetscScalar *bl; 

    ot::DA* da = m_ts->getMass()->getOctDA();

    unsigned int maxD = da->getMaxDepth();
    unsigned int balOctmaxD = maxD - 1;

    for (unsigned int g=0; g<m_radialBasis.size(); g++) {
      for (int b=0; b<knotsize; b++) {
        da->vecGetBuffer(ibldt[b], bl, false, true, false, 3);
        for ( da->init<ot::DA::ALL>(), da->init<ot::DA::WRITABLE>(); da->curr() < da->end<ot::DA::ALL>(); da->next<ot::DA::ALL>()) {
          unsigned int i = da->curr();
          Point pt;
          pt = da->getCurrentOffset();
          unsigned levelhere = da->getLevel(da->curr()) - 1;
          double hxOct = (double)((double)(1<<(balOctmaxD - levelhere))/(double)(1 << balOctmaxD));

          double x = (double)(pt.xint())/((double)(1<<(maxD-1)));
          double y = (double)(pt.yint())/((double)(1<<(maxD-1)));
          double z = (double)(pt.zint())/((double)(1<<(maxD-1)));

          Point px(x,y,z);

          gx = m_radialBasis[g].getValue(px);

          sendVec[g*bsize + 3*b]     += bl[3*i]*gx;
          sendVec[g*bsize + 3*b + 1] += bl[3*i+1]*gx;
          sendVec[g*bsize + 3*b + 2] += bl[3*i+2]*gx;
        }
        da->vecRestoreBuffer(ibldt[b], bl, false, true, false, 3);
      } // b
    } // g
  }

  delete [] currBasis;
  MPI_Barrier(MPI_COMM_WORLD);

  // update copies on all ...
  MPI_Allreduce ( sendVec, pVec, sz,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

  VecRestoreArray(params, &pVec);

  delete [] sendVec;
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}

#endif
