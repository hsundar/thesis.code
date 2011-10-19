#ifndef __INVERSE_SOLVER_H_
#define __INVERSE_SOLVER_H_

#include "petscksp.h"
#include "petscda.h"
#include "petscdmmg.h"
#include "timeStepper.h"
#include "stsdamgHeader.h"
// #include "rpHeader.h"


class inverseSolver {

  public:
    // Constructor and Destructors
    inverseSolver();
    //  inverseSolver(timeStepper* ts1);
    virtual ~inverseSolver();

    // Operations

    /** @name access methods **/ 
    //@{
    // set Parameters

    void setRegularizationParameter(double beta)
    {
      m_beta = beta;
    }
    void setParamsTolerance (double x)
    {
      m_paramsTolerance = x;
    } 
    void setFunctionTolerance (double x)
    {
      m_functionTolerance = x;
    } 
    double getFunctionTolerance () const
    {
      return m_functionTolerance;
    }
    int getMaximumNumberOfIterations () const
    {
      return m_maxIterations;
    }
    int getFinalNumberOfIterations () const
    {
      return m_numIterations;
    }
    double getCostFunctionValue () const
    {
      return m_costFunctionValue;
    }

    // Partial Observations.
    PetscErrorCode setPartialObservations(Vec pObs) {
      m_vecPartialObservations = pObs;
      m_bUsePartialObservations = true;
      return(0);
    }

    // 
    PetscErrorCode setInitialGuess(Vec initialGuess);

    // 
    PetscErrorCode setTimeStepper(timeStepper* ts);

    // set right hand sides
    PetscErrorCode setForwardRhs(feVec* forwardRhs);
    PetscErrorCode setAdjointRhs(Vec adjointRhs);

    //  PetscErrorCode setReducedGradient(feVec* reducedGradient);
    //  PetscErrorCode setReducedHessian(feMat* reducedHessian);

    // Get control, adjoint and State
    PetscErrorCode getControlStep(Vec controlStep);
    PetscErrorCode getCurrentControl(Vec currentControl);

    PetscErrorCode getCurrentState(Vec currentState);
    PetscErrorCode getCurrentAdjoint(Vec currentAdjoint);

    bool isOptimizing ()const
    {
      return m_isOptimizing;
    };

    virtual void  hessianMatMult(Vec _in, Vec _out)= 0;

    virtual void mghessianMatMult(DA _da, Vec _in, Vec _out) = 0;

    static PetscErrorCode MatMult(Mat M, Vec In, Vec Out){
#ifdef __DEBUG__    
      std::cout << "Entering " << __func__ << std::endl;
#endif    
      inverseSolver *contxt;
      MatShellGetContext(M,(void**)&contxt);

      contxt->hessianMatMult(In,Out);
      return(0);
    }

    static PetscErrorCode MGMatMult(Mat M, Vec In, Vec Out){
      stsDMMG *contxt;
      MatShellGetContext(M,(void**)&contxt);

      DA da = (DA)(((stsDMMG)contxt)->dm);

      ((inverseSolver*)(((stsDMMG)contxt)->user))->mghessianMatMult(da,In,Out);

#ifdef __DEBUG__
      PetscPrintf(0,"Exiting MGMatMult\n");
#endif
      return(0);
    }

    static PetscErrorCode CreateHessian(stsDMMG dmmg, Mat *J){
      DA da = (DA)(dmmg->dm);
      int m,n,xm,ym,zm;
      int ierr;
      DALocalInfo info;

      ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
      ierr = DAGetCorners(da,0,0,0,&xm,&ym,&zm); CHKERRQ(ierr);
      m=n=xm*ym*zm;

      std::cout << "size @ level = "<< m << std::endl;
      ierr = MatCreateShell(PETSC_COMM_WORLD,m*info.dof,n*info.dof,PETSC_DETERMINE,PETSC_DETERMINE,dmmg,J); CHKERRQ(ierr);
      ierr = MatShellSetOperation(*J,MATOP_MULT,(void(*)(void))MGMatMult); CHKERRQ(ierr);

      return(0);
    }

    static PetscErrorCode ComputeHessian(stsDMMG dmmg, Mat A, Mat B){

      return(0);
    }

    static PetscErrorCode ComputeRHS(stsDMMG dmmg, Vec b){

      // int ierr;
      // int size;

      VecCopy(((inverseSolver*)dmmg->user)->m_vecReducedGradient,b);
      return(0);
    }

    //@}
  protected:
    /** @name Data members **/ 
    //@{
    double m_costFunctionValue;
    double m_paramsTolerance;
    double m_functionTolerance;
    double m_beta;

    int m_maxIterations;
    int m_numIterations;
    bool m_isOptimizing;

    // Timestepper
    timeStepper* m_ts;

    // Initial Guess
    Vec m_vecInitialControl;

    // For dealing for partial observations
    Vec m_vecPartialObservations;
    bool m_bUsePartialObservations;

    // Working Vectors
    Vec m_vecCurrentControl;
    Vec m_vecCurrentState;
    Vec m_vecCurrentAdjoint;
    // control Step 
    Vec m_vecControlStep;

    // Right hand side objects
    feVec* m_vecForwardRhs;
    Vec m_vecAdjointRhs;

    // Reduced Gradient
    Vec m_vecReducedGradient;

    // Reduced Hessian
    Mat m_matReducedHessian;

    // Krylov solver for the step
    KSP m_ksp;

    // Multigrid solver for the step
    stsDMMG *m_dmmg;

};

#endif 
