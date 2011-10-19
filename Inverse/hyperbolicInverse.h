/**
 * @file hyperbolicInverse.h
 * @brief Main class for an inverse hyperbolic problem
 * @author Hari Sundar
 * @date   8/13/07
 * 
 * Main class for inverse hyperbolic problem.
 * 
 **/

#ifndef _HYPERBOLIC_INVERSE_H_
#define _HYPERBOLIC_INVERSE_H_

#include "inverseSolver.h"
#include "femUtils.h"

class hyperbolicInverse : public inverseSolver { // hyperbolic inverse

public:

  hyperbolicInverse() {}

  virtual ~hyperbolicInverse() {}

  virtual int destroy() {
    // Allocate memory for working vectors
    CHKERRQ(VecDestroy(m_vecCurrentControl));
    CHKERRQ(VecDestroy(m_vecControlStep));
    CHKERRQ(VecDestroy(m_vecReducedGradient));

    CHKERRQ(MatDestroy(m_matReducedHessian));
    CHKERRQ(KSPDestroy(m_ksp));
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
  virtual void mghessianMatMult(DA _da, Vec _in, Vec _out) {}

  PetscErrorCode setObservations(std::vector<Vec> obs) {
    m_vecObservations = obs;
    return(0);
  }

protected:
  // In case of time dependent problems, it's faster to save the timesteps separately.
  std::vector<Vec> m_vecObservations;

  Vec m_vecForwardInitialDisplacement;
  Vec m_vecForwardInitialVelocity;

};

/**
 *	@brief The initialization function where the reduced Hessian matrix, KSP context are created
 * @return 0 if successful, 1 otherwise
 * Matrix of shell type is created which does only a matvec
 **/
int hyperbolicInverse::init()
{
  // std::cout << RED"Entering "NRM << __func__ << std::endl;
  int matsize;
  int ierr;

  // Get reduced Hessian matrix size
  ierr = VecDuplicate(m_vecInitialControl,&m_vecCurrentControl); CHKERRQ(ierr);
  ierr = VecDuplicate(m_vecInitialControl,&m_vecReducedGradient); CHKERRQ(ierr);
  ierr = VecDuplicate(m_vecInitialControl,&m_vecControlStep); CHKERRQ(ierr);
  ierr = VecGetLocalSize(m_vecInitialControl,&matsize); CHKERRQ(ierr);

  // Create Reduced Hessian
  ierr = MatCreateShell(PETSC_COMM_WORLD,matsize,matsize,PETSC_DETERMINE,PETSC_DETERMINE,this,&m_matReducedHessian); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_matReducedHessian, MATOP_MULT,(void(*)(void))(MatMult)); CHKERRQ(ierr);

  // Create a KSP context to solve for the reduced Hessian system
  ierr = KSPCreate(PETSC_COMM_WORLD,&m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(m_ksp,"inv_"); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_ksp, m_matReducedHessian, m_matReducedHessian, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetType(m_ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);

  return(0);
}

int hyperbolicInverse::solve()
{
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif
  
  int ierr;
  double norm;

  // Set the initial guess to the current control
  ierr = VecCopy(m_vecInitialControl, m_vecCurrentControl); CHKERRQ(ierr);

#ifdef __DEBUG__
  VecNorm(m_vecCurrentControl,NORM_2,&norm);
  PetscPrintf(0,"norm of initial guess = %g\n",norm);
#endif

  // initiate the step to zero
  ierr = VecZeroEntries(m_vecControlStep); CHKERRQ(ierr);

  // compute the reduced Gradient
  PetscPrintf(0,"Setting Reduced Gradient\n");
  setReducedGradient();

  // Solve for the step using the reduced Hessian
  PetscPrintf(0,"Hessian KSP Solve\n");
  ierr = KSPSolve(m_ksp, m_vecReducedGradient, m_vecControlStep); CHKERRQ(ierr);
  PetscPrintf(0,"Done Hessian KSP Solve\n");

  PetscReal rnorm;
  PetscInt its;
  KSPGetResidualNorm(m_ksp, &rnorm);
  KSPGetIterationNumber(m_ksp,&its);
  // Print the final residual ...
  PetscPrintf(0, "Final residual norm is %g\n", rnorm);
  PetscPrintf(0, "Total number of iterations: %d\n", its);

  /*
  Mat Hessian, HTrans;
  double hnorm;
  ierr = KSPComputeExplicitOperator(m_ksp,&Hessian); CHKERRQ(ierr);
  ierr = MatTranspose(Hessian, &HTrans);
  MatAXPY(Hessian, -1.0, HTrans, DIFFERENT_NONZERO_PATTERN);
  MatNorm(Hessian, NORM_INFINITY, &hnorm);

  PetscPrintf(0," HNorm = %g\n", hnorm);
  */
  // new control = old control - step
  ierr = VecAXPY(m_vecCurrentControl, -1.0, m_vecControlStep); CHKERRQ(ierr);

#ifdef __DEBUG__
  VecNorm(m_vecCurrentControl,NORM_INFINITY,&norm);
  PetscPrintf(0,"norm of solution = %g\n",norm);
#endif
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
  return(0);
}

/**
 *	@brief set the reduced Gradient for the current problem, this function does a forward solve and then -C^T J^-T( y - y*)
 * 
 * Two full solves (one forward + one adjoint)
 **/
bool hyperbolicInverse::setReducedGradient() {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif

  double norm;

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

  timeInfo *ti = ts->getTimeInfo();
  unsigned NT = (int)(ceil((ti->stop - ti->start)/ti->step));

  // Get the force from the timestepper
  Fdynamic = (cardiacDynamic *)ts->getForce();

  // Set the force to the current control
  std::vector<Vec> currControl;
  splitVec(m_vecCurrentControl, currControl, NT+1);
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

  /*
  concatenateVecs(solvec, m_vecReducedGradient, false);
  VecNorm(m_vecReducedGradient, NORM_2, &norm);
  PetscPrintf(0, "Solution Norm is %g\n", norm);
  */
  // VecGetLocalSize( m_vecReducedGradient, &lSize );
  // std::cout << CYN"After fwd: Reduced gradient size is "NRM << lSize << std::endl;

  // std::cout << "Sizes are " << solvec.size() << ", " << m_vecObservations.size() << std::endl;

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

  // clear memory ...
  for (int i=0; i<m_vecObservations.size(); i++) {
    if (m_vecObservations[i] != NULL) {
      VecDestroy(m_vecObservations[i]);
    }
  }

  // concatenateVecs(m_vecAdjoints, m_vecReducedGradient, false);
  // PetscScalar norm;
  // VecNorm(m_vecReducedGradient, NORM_2, &norm);
  // PetscPrintf(0, "Adjoints Norm is %g\n", norm);
  // concatenateVecs(solvec, m_vecReducedGradient, false);
  // VecNorm(m_vecReducedGradient, NORM_2, &norm);
  // PetscPrintf(0, "Final RHS Norm is %g\n", norm);

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
  for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }

  // get the solution
  solvec = ts->getSolution();
  //VecCopy(solvec[0],m_vecReducedGradient);

  // Now concatenate the solution ...
  concatenateVecs(solvec, m_vecReducedGradient, false);

  // clear memory ...
  for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }

  // scale the reduced Gradient
  VecScale(m_vecReducedGradient,-1.0);
  // VecGetLocalSize( m_vecReducedGradient, &lSize );
  // std::cout << YLW"after Adj: Reduced gradient size is "NRM << lSize << std::endl;

  // add the contribution of the regularization parameter
  VecAXPY(m_vecReducedGradient, m_beta, m_vecCurrentControl);

//#ifdef __DEBUG__
  // VecNorm(m_vecReducedGradient,NORM_2,&norm);
  // PetscPrintf(0,"norm of reduced Gradient = %f\n", norm);

  VecDestroy(initD);
  VecDestroy(initV);

// #endif
  // std::cout << GRN"Leaving "NRM << __func__ << std::endl;

  return true;
}
/**
 *	@brief set the reduced Hessian matrix vector product, this function does one forward solve and one adjoint solve
 *
 * Two full solves : one forward + one adjoint
 **/
void hyperbolicInverse::hessianMatMult(Vec In, Vec Out)
{
  // double norm;
  // std::cout << RED"Entering "NRM << __func__ << std::endl;
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

  timeInfo *ti = ts->getTimeInfo();
  unsigned NT = (int)(ceil((ti->stop - ti->start)/ti->step));

  // Get the force from the timestepper
  Fdynamic = (cardiacDynamic *)ts->getForce();

  // Set the force to In
  std::vector<Vec> currControl;
  splitVec(In, currControl, NT+1);
  Fdynamic->setFDynamic(currControl);

  ts->setForceVector(Fdynamic);
  ts->clearMonitor();
  // ts->init();
  ts->solve();

  // state variable
  std::vector<Vec> solvec;
  solvec = ts->getSolution();

  Vec tmp;
  VecDuplicate(m_vecReducedGradient, &tmp);
  concatenateVecs(solvec, tmp, false);
  // VecNorm(tmp, NORM_2, &norm);
  // PetscPrintf(0, "Sol Norm after HessFwd is %g\n", norm);

  // Scale to the set the right hand side of adjoint
  for (unsigned int i=0; i<solvec.size(); i++) {
      VecScale(solvec[i],-1.0);
      if (m_bUsePartialObservations) {
        VecPointwiseMult(solvec[i], solvec[i], m_vecPartialObservations);
      }
  }

  // concatenateVecs(solvec, tmp, false);
  //VecNorm(tmp, NORM_2, &norm);
  // PetscPrintf(0, "Sol Norm after HessFwd Scale is %g\n", norm);


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

  // ts->init();
  // adjoint solve
  ts->solve();

  // clear memory ...
  for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }

  // Adjoint variable
  solvec = ts->getSolution();

  // Now concatenate the solution ...
  concatenateVecs(solvec, Out, false);

  // clear memory ...
  for (int i=0; i<solvec.size(); i++) {
    if (solvec[i] != NULL) {
      VecDestroy(solvec[i]);
    }
  }

  // scaling the adjoint variable
  VecScale(Out, -1.0);

  //VecNorm(Out, NORM_2, &norm);
  //PetscPrintf(0, "Norm after HessAdj is %g\n", norm);

  // add the contribution of the regularization parameter
  VecAXPY(Out, m_beta, In);
  // std::cout << GRN"Leaving "NRM << __func__ << std::endl;
  VecDestroy(initD);
  VecDestroy(initV);
}

#endif
