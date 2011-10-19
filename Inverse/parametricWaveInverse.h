/**
 * @file parametricWaveInverse.h
 * @brief Main class for an inverse scalar Hyperbolic problem
 * @author Hari Sundar 
 * @date   11/13/07 
 * 
 * Main class for inverse scalar Hyperbolic problem where the 
 * force is parametrised. 
 * 
 **/

#ifndef _PARAMETRIC_WAVE_INVERSE_H_
#define _PARAMETRIC_WAVE_INVERSE_H_

#include <vector> 
#include "inverseSolver.h"
#include "femUtils.h"

class parametricWaveInverse : public inverseSolver {

public:

  parametricWaveInverse() { m_bComputeBasisOnTheFly = true;}

  ~parametricWaveInverse() {}


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

  PetscErrorCode setAdjoints(std::vector<Vec> adjoints) {
    m_vecAdjoints = adjoints;
    return(0);
  }

  PetscErrorCode setForceBasis(std::vector<std::vector< Vec > > basis) {
    m_forceBasis = basis;
    m_bComputeBasisOnTheFly = false;
    return(0);
  }

  PetscErrorCode setNumberOfParameter(unsigned int numP) {
    m_iNumParams = numP;
    m_bComputeBasisOnTheFly = true;
    return(0);
  }

  // Functions to hanlde the full / parametrized representations
  void getForces(Vec params, std::vector<Vec> &forces);
  void getParams(std::vector<Vec> forces, Vec params);


protected:
  // In case of time dependent problems, it's faster to save the timesteps separately.
  std::vector<Vec> m_vecAdjoints;

  // The parametrized force vectors ...
  std::vector< std::vector <Vec> > m_forceBasis;
  unsigned int m_iNumParams;

  bool m_bComputeBasisOnTheFly;

  Vec m_vecForwardInitialDisplacement;
  Vec m_vecForwardInitialVelocity;

};

/**
 *	@brief The initialization function where the reduced Hessian matrix, KSP context are created
 * @return 0 if successful, 1 otherwise
 * Matrix of shell type is created which does only a matvec
 **/
int parametricWaveInverse::init()
{
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

  // Create Reduced Hessian
  // ierr = MatCreateSeqDense(PETSC_COMM_SELF, matsize, matsize, PETSC_NULL, &m_matReducedHessian); CHKERRQ(ierr);
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

int parametricWaveInverse::solve()
{

  // std::cout << RED"Entering "NRM << __func__ << std::endl;
  int ierr;
  // double norm;

  // Set the initial guess to the current control
  ierr = VecCopy(m_vecInitialControl, m_vecCurrentControl); CHKERRQ(ierr);

#ifdef __DEBUG__
  VecNorm(m_vecCurrentControl,NORM_2,&norm);
  PetscPrintf(0,"norm of initial guess = %g\n",norm);
#endif

  // initiate the step to zero
  ierr = VecZeroEntries(m_vecControlStep); CHKERRQ(ierr);

  // compute the reduced Gradient
  // PetscPrintf(0,"Setting Reduced Gradient\n");

  setReducedGradient();

  /* 
   PetscInt lSize=0;
   VecGetLocalSize( m_vecControlStep, &lSize );
   std::cout << BLU"Control step size is "NRM << lSize << std::endl;
   VecGetLocalSize( m_vecReducedGradient, &lSize );
   std::cout << BLU"Reduced gradient size is "NRM << lSize << std::endl;
 */

  // Solve for the step using the reduced Hessian
  // PetscPrintf(0,"Hessian KSP Solve\n");
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
  // std::cout << GRN"Leaving "NRM << __func__ << std::endl;

  return(0);
}

/**
 *	@brief set the reduced Gradient for the current problem, this function does a forward solve and then -C^T J^-T( y - y*)
 * 
 * Two full solves (one forward + one adjoint)
 **/
bool parametricWaveInverse::setReducedGradient()
{
  // std::cout << RED"Entering "NRM << __func__ << std::endl;
  // double norm;

  // Initiate the reduced Gradient to zero
  VecZeroEntries(m_vecReducedGradient);

  // Fdynamic which is the control 
  fdynamicVector *Fdynamic = new fdynamicVector(feVec::PETSC);

  // Solver 
  newmark *ts = new newmark;

  ts->setInitialDisplacement(m_vecForwardInitialDisplacement);
  ts->setInitialVelocity(m_vecForwardInitialVelocity);

  // Forward solve
  ts = (newmark *)m_ts;
  ts->setAdjoint(false);

  ts->setTimeFrames(1);
  ts->storeVec(true);

  // timeInfo *ti = ts->getTimeInfo();
  // unsigned NT = (int)(ceil((ti->stop - ti->start)/ti->step));

  // Get the force from the timestepper
  Fdynamic = (fdynamicVector *)ts->getForce();

  // Set the force to the current control
  std::vector<Vec> currControl;
  // splitVec(m_vecCurrentControl, currControl, NT+1);
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

   // Now can clear memory ...
  for (int i=0; i<currControl.size(); i++) {
    if (currControl[i] != NULL) {
      VecDestroy(currControl[i]);
    }
  }

  /*
  concatenateVecs(solvec, m_vecReducedGradient, false);
  VecNorm(m_vecReducedGradient, NORM_2, &norm);
  PetscPrintf(0, "Solution Norm is %g\n", norm);
  */

  // VecGetLocalSize( m_vecReducedGradient, &lSize );
  // std::cout << CYN"After fwd: Reduced gradient size is "NRM << lSize << std::endl;

#ifdef __DEBUG__
  PetscPrintf(0,"size of solvec is %d\n", solvec.size());
  VecNorm(solvec[0],NORM_INFINITY,&norm);
  PetscPrintf(0,"norm of state in reduced gradient = %f\n", norm);
#endif

  // std::cout << solvec.size() << ", " << m_vecAdjoints.size() << std::endl;

  // Set the right hand side of the adjoint, (state - data)
  for (unsigned int i=0; i<solvec.size(); i++) {
    VecAYPX(solvec[i], -1.0, m_vecAdjoints[i]);  
    if (m_bUsePartialObservations) {
      VecPointwiseMult(solvec[i], solvec[i], m_vecPartialObservations);
    }
  }

  // Now can clear memory ...
  for (int i=0; i<m_vecAdjoints.size(); i++) {
    if (m_vecAdjoints[i] != NULL) {
      VecDestroy(m_vecAdjoints[i]);
    }
  }

  /*
  concatenateVecs(m_vecAdjoints, m_vecReducedGradient, false);
  //PetscScalar norm;
  VecNorm(m_vecReducedGradient, NORM_2, &norm);
  PetscPrintf(0, "Adjoints Norm is %g\n", norm);
  concatenateVecs(solvec, m_vecReducedGradient, false);
  VecNorm(m_vecReducedGradient, NORM_2, &norm);
  PetscPrintf(0, "Final RHS Norm is %g\n", norm);
  */

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
  // concatenateVecs(solvec, m_vecReducedGradient, false);

  getParams(solvec, m_vecReducedGradient);

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

#ifdef __DEBUG__
  VecNorm(m_vecReducedGradient,NORM_2,&norm);
  PetscPrintf(0,"norm of reduced Gradient = %f\n", norm);
#endif
  VecDestroy(initD);
  VecDestroy(initV);

  // std::cout << GRN"Leaving "NRM << __func__ << std::endl;

  return true;
}
/**
 *	@brief set the reduced Hessian matrix vector product, this function does one forward solve and one adjoint solve
 *
 * Two full solves : one forward + one adjoint
 **/
void parametricWaveInverse::hessianMatMult(Vec In, Vec Out)
{
  // double norm;
  // std::cout << RED"Entering "NRM << __func__ << std::endl;
  VecZeroEntries(Out);

  newmark *ts = new newmark;

  // Fdynamic which is the control 
  fdynamicVector *Fdynamic = new fdynamicVector(feVec::PETSC);

  // Forward solve
  ts = (newmark *)m_ts;

  ts->setInitialDisplacement(m_vecForwardInitialDisplacement);
  ts->setInitialVelocity(m_vecForwardInitialVelocity);

  ts->setTimeFrames(1);
  ts->storeVec(true);
  ts->setAdjoint(false);

  // timeInfo *ti = ts->getTimeInfo();
  // unsigned NT = (int)(ceil((ti->stop - ti->start)/ti->step));

  // Get the force from the timestepper
  Fdynamic = (fdynamicVector *)ts->getForce();

  // Set the force to In
  std::vector<Vec> currControl;
  // splitVec(In, currControl, NT+1);
  getForces(In, currControl);
  Fdynamic->setFDynamic(currControl);

  ts->setForceVector(Fdynamic);
  ts->clearMonitor();
  // ts->init();
  ts->solve();

  // state variable
  std::vector<Vec> solvec;
  solvec = ts->getSolution();

   // Now can clear memory ...
  for (int i=0; i<currControl.size(); i++) {
    if (currControl[i] != NULL) {
      VecDestroy(currControl[i]);
    }
  }

  // Vec tmp;
  // VecDuplicate(m_vecReducedGradient, &tmp);
  // concatenateVecs(solvec, tmp, false);
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
  // VecNorm(tmp, NORM_2, &norm);
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
  // concatenateVecs(solvec, Out, false);
  getParams(solvec, Out);

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

// Functions to hanlde the full / parametrized representations
void parametricWaveInverse::getForces(Vec params, std::vector<Vec> &forces) {
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

  if (!m_bComputeBasisOnTheFly) {
    // create and initialize to alpha_0 * Basis_0
    for (unsigned int i=0; i<m_forceBasis[0].size(); i++) {
      Vec tmp;
      VecDuplicate(m_forceBasis[0][i], &tmp);
      VecZeroEntries(tmp);
      VecAXPY(tmp, pVec[0], m_forceBasis[0][i]);
      forces.push_back(tmp);
    }

    for (unsigned int i=1; i<m_forceBasis.size(); i++) {
      for (unsigned int j=0; j<m_forceBasis[i].size(); j++) {
        VecAXPY(forces[j], pVec[i], m_forceBasis[i][j]);
      }
    }
  } else {
    DA da = m_ts->getMass()->getDA();
    timeInfo *ti = m_ts->getTimeInfo();


    unsigned int numSteps = (unsigned int)(ceil(( ti->stop - ti->start)/ti->step));
    // create and initialize to 0
    for (unsigned int i=0; i<numSteps+1; i++) {
      Vec tmp;
      DACreateGlobalVector(da, &tmp);
      VecZeroEntries(tmp);
      forces.push_back(tmp);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PetscScalar ***tauArray;

    unsigned int paramTimeSteps = (unsigned int)(ceil(( (double)(numSteps))/ ((double)(2*m_iNumParams)) ));
    double acx,acy,acz;

    int x, y, z, m, n, p;
    int mx,my,mz;

    DAGetCorners(da, &x, &y, &z, &m, &n, &p);
    DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0);

    double hx = 1.0/(mx-1.0);

    for (int b=0; b<m_iNumParams; b++) {
      std::vector<Vec> tau;
      unsigned int tBegin = paramTimeSteps*b;
      unsigned int tEnd   = tBegin + numSteps/2; // paramTimeSteps*(b+2);

      // std::cout << "For param " << b << ": Time step range is " << tBegin << " -> " << tEnd << std::endl; 
      for (unsigned int t=0; t<numSteps+1; t++) {
        double newTime = (ti->step*(t-tBegin)*numSteps)/((double)(paramTimeSteps));

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
  }

  VecRestoreArray(params, &pVec);
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}

void parametricWaveInverse::getParams(std::vector<Vec> forces, Vec params) {
#ifdef __DEBUG__
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif

  // To recover the parameters, given the force, we need to integrate over space and time
  // alpha_i = \int \int F B_i

  // Need the DA for this ...
  DA da = m_ts->getMass()->getDA();
  timeInfo *ti = m_ts->getTimeInfo();

  unsigned int numSteps = (unsigned int)(ceil(( ti->stop - ti->start)/ti->step));

  // Get spatial dimentions and the like ...
  int x, y, z, m, n, p;
  int mx,my,mz;
  int sz;

  DAGetCorners(da, &x, &y, &z, &m, &n, &p); 
  DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,0,0,0,0); 

  double hx = 1.0/((double)mx-1.0);
  double dt = ti->step;

  VecZeroEntries(params);
  PetscScalar * pVec;
  VecGetArray(params, &pVec);
  VecGetLocalSize(params, &sz);
  PetscScalar *sendVec = new PetscScalar[sz];

  for (int i=0; i<sz; i++) {
    sendVec[i] = 0.0;
  }

  double scaleT, scaleX, scaleY, scaleZ;
  PetscScalar ***basis, ***lambda;

  unsigned int numP=0;

  if (m_bComputeBasisOnTheFly) {
    numP = m_iNumParams;
  } else {
    numP = m_forceBasis.size();
  }

  unsigned int paramTimeSteps = (unsigned int)(ceil(( (double)(numSteps))/ ((double)(2*m_iNumParams)) ));
  double acx,acy,acz;

  for (unsigned int b=0; b<numP; b++) {
    unsigned int tBegin = paramTimeSteps*b;
    unsigned int tEnd   = tBegin + numSteps/2;
    for (unsigned int t=0; t<numSteps+1; t++) {
      if (!m_bComputeBasisOnTheFly) {
        DAVecGetArray(da, m_forceBasis[b][t], &basis);
      }

      DAVecGetArray(da, forces[t], &lambda);

      if (!t || (t == numSteps))
        scaleT = 0.5;
      else
        scaleT = 1.0;
      for (int k=z; k<z+p; k++) {
        if (!k || (k == (mz-1)))
          scaleZ = 0.5;
        else
          scaleZ = 1.0;
        for (int j=y; j<y+n; j++) {
          if (!j || (j == (my-1)))
            scaleY = 0.5;
          else
            scaleY = 1.0;
          for (int i=x; i<x+m; i++) {
            if (!i || (i == (mx-1)))
              scaleX = 0.5;
            else
              scaleX = 1.0;
            // The actual product ...
            if (!m_bComputeBasisOnTheFly) {
              sendVec[b] += scaleT*scaleZ*scaleY*scaleX*basis[k][j][i]*lambda[k][j][i];
            } else {
              if ( (t>=tBegin) && (t<=tEnd)) {
                double newTime = (dt*(t-tBegin)*numSteps)/((double)(paramTimeSteps));
                acx = (i)*hx; acy = (j)*hx; acz = (k)*hx;
                sendVec[b] += scaleT*scaleZ*scaleY*scaleX*lambda[k][j][i]*sin(M_PI*newTime)*cos(2*M_PI*acx)*cos(2*M_PI*acy)*cos(2*M_PI*acz);
              }
            }
          }
        }
      }
      if (!m_bComputeBasisOnTheFly) {
        DAVecRestoreArray(da, m_forceBasis[b][t], &basis);
      }
      DAVecRestoreArray(da, forces[t], &lambda);
    }
    sendVec[b] *= hx*hx*hx*dt;
  }

  // update copies on all ...
  MPI_Allreduce ( sendVec, pVec, sz,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

  VecRestoreArray(params, &pVec);
  delete [] sendVec;
  // VecView(params, 0);
#ifdef __DEBUG__
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif
}

#endif
