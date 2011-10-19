#include "inverseSolver.h"

// constructor
inverseSolver::inverseSolver()
{
  m_bUsePartialObservations = false;
}
// destructor
inverseSolver::~inverseSolver()
{
  
}

// set Timestepper
PetscErrorCode inverseSolver::setTimeStepper(timeStepper *ts)
{
  m_ts = ts;
  return(0);
}

// set InitialGuess

PetscErrorCode inverseSolver::setInitialGuess(Vec InitialGuess)
{
  int ierr;
  ierr = VecDuplicate(InitialGuess,&m_vecInitialControl); CHKERRQ(ierr);
  ierr = VecCopy(InitialGuess,m_vecInitialControl); CHKERRQ(ierr);
  return(0);
}

// set Forward rhs
PetscErrorCode inverseSolver::setForwardRhs(feVec* forwardRhs)
{
  m_vecForwardRhs = forwardRhs;
  return(0);
}

// set Adjoint rhs
PetscErrorCode inverseSolver::setAdjointRhs(Vec adjointRhs)
{
  m_vecAdjointRhs = adjointRhs;
  return(0);
}

// get control step Copy control
PetscErrorCode inverseSolver::getControlStep(Vec controlStep)
{
  int ierr;
  ierr = VecCopy(m_vecControlStep,controlStep); CHKERRQ(ierr);

  return(0);
}

// get currentState
PetscErrorCode inverseSolver::getCurrentState(Vec currentState)
{
  int ierr;
  ierr = VecCopy(m_vecCurrentState,currentState); CHKERRQ(ierr);

  return(0);
}

// get adjoint variable
PetscErrorCode inverseSolver::getCurrentAdjoint(Vec currentAdjoint)
{
  int ierr;
  ierr = VecCopy(m_vecCurrentAdjoint, currentAdjoint); CHKERRQ(ierr);

  return(0);
}

// get adjoint variable
PetscErrorCode inverseSolver::getCurrentControl(Vec currentControl)
{
  int ierr;
  ierr = VecCopy(m_vecCurrentControl, currentControl); CHKERRQ(ierr);

  return(0);
}

// 
