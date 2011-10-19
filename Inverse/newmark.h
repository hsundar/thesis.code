/**
*  @file	Newmark.h
*  @brief	Implements the Newmark scheme for solving the dynamic hyperbolic problem.
*  @author	Hari Sundar
*  @date	5/12/07
*
*  Implements the Newmark scheme for solving the dynamic hyperbolic problem. The problems
*  are assumed to be of the type,
* 
*  \f[
*  		\bf M\ddot{d} + C\dot{d} + Kd = F
*  \f]
* 
*  It is important to specify the mass, damping and stiffness matrices and a force vector. The matrices
*  need to be derived from feMatrix and the force vector should be derived from feVector. The 
*  time step and the duration of the problem also need to be specified.
* 
*  The problem needs the initial displacement \f${\bf u}_0\f$ and the initial velocity, \f${\bf v}_0\f$
*  to be specified. The final displacement fields, i.e., the motion estimates can be obtained once
*  the solution has converged.   
**/

#ifndef _NEWMARK_H_
#define _NEWMARK_H_

#include <fstream>
#include <iomanip>

#include "timeStepper.h"
#include "colors.h"


/**
*  @brief	Implements the Newmark scheme for solving the dynamic hyperbolic problem.
*  @author	Hari Sundar
*  @date	5/12/07
*
*  Implements the Newmark scheme for solving the dynamic hyperbolic problem. The problems
*  are assumed to be of the type,
* 
*  \f[
*  		\bf M\ddot{d} + C\dot{d} + Kd = F
*  \f]
* 
*  It is important to specify the mass, damping and stiffness matrices and a force vector. The matrices
*  need to be derived from feMatrix and the force vector should be derived from feVector. The 
*  time step and the duration of the problem also need to be specified.
* 
*  The problem needs the initial displacement \f${\bf u}_0\f$ and the initial velocity, \f${\bf v}_0\f$
*  to be specified. The final displacement fields, i.e., the motion estimates can be obtained once
*  the solution has converged.   
**/

class newmark : public timeStepper {
public:


	virtual int init();
	virtual int solve();

	virtual int destroy() {
		// Allocate memory for working vectors
		CHKERRQ(VecDestroy(m_vecSolution));
		CHKERRQ(VecDestroy(m_vecVelocity));
		CHKERRQ(VecDestroy(m_vecAccn));

		CHKERRQ(VecDestroy(m_vecRHS));

		CHKERRQ(VecDestroy(m_vec_du));
		CHKERRQ(VecDestroy(m_vec_dv));
		CHKERRQ(VecDestroy(m_vec_da));

		CHKERRQ(MatDestroy(m_matJacobian));
		CHKERRQ(KSPDestroy(m_ksp));

		CHKERRQ(MatDestroy(m_matMass));
		CHKERRQ(KSPDestroy(m_AccnKSP));

		return true;
	}

	/**
	 *	@brief The Jacobian matmult routine by matrix free 
	 *	       method using the given stiffness, damping and mass matrices
	 *  @param m_matJacobian PETSC Matrix, this is of type matrix shell
	 *  @param In  PETSC Vec, input vector
	 *  @param Out PETSC Vec, Output vector
	 **/
	virtual void jacobianMatMult(Vec In, Vec Out);
	virtual void massMatMult(Vec In, Vec Out);

	virtual void jacobianGetDiagonal(Vec diag);
	virtual void massGetDiagonal(Vec diag);

	/**
	 *	@brief Sets the right hand side vector given the current solution
	 *  @param m_vecCurrentSolution PETSC Vec, current solution
	 *  @param m_vecNextRHS PETSC Vec, next right hand side
	 **/
	virtual bool setRHS();

	// Newmark Solve for initial accleration

	static void InitMatMult(Mat M, Vec In, Vec Out) {
		newmark *contxt;
		MatShellGetContext(M, (void**)&contxt);

		contxt->massMatMult(In,Out);
	}
	static void InitMatGetDiagonal(Mat M, Vec diag) {
		newmark *contxt;
		MatShellGetContext(M, (void**)&contxt);

		contxt->massGetDiagonal(diag);
	}

	virtual bool setAccnRHS();


	// Not yet implemented ...
	void  mgjacobianMatMult(DA _da, Vec _in, Vec _out) {
	};

	bool  setRHSFunction(Vec _in, Vec _out) {
		return true;
	};

	/**
	*  @brief set time frames at which solution is stored
	**/
	virtual int setTimeFrames(int mon) {
		m_iMon = mon;
		return 0;
	}

	/**
	* @brief observer operator which saves the solution after every t timesteps
	**/
	virtual int monitor();
	int clearMonitor() {
		/*
		for (int i=0; i<m_solVector.size(); i++) {
			if (m_solVector[i] != NULL) {
				VecDestroy(m_solVector[i]);
			}
		} */
		m_solVector.clear();
		CHKERRQ( VecZeroEntries( m_vecSolution ) );
		CHKERRQ( VecZeroEntries( m_vecVelocity ) );
		CHKERRQ( VecZeroEntries( m_vecAccn ) );
		return(0);
	}

	/**
	* @brief get the solution vector
	**/
	virtual std::vector<Vec> getSolution() {
		return m_solVector;
	}

	void storeVec(bool flag) {
		m_bStoreVec = flag;
	}

	void useMatrixFree(bool mfree) {
		m_bMatrixFree = mfree;
	}

	void damp(bool f) {
		m_bDamp = f;
	}

protected:
	double    m_dBeta;
	double    m_dGamma;

	bool      m_bDamp;

	Vec       m_vec_du;
	Vec       m_vec_dv;
	Vec       m_vec_da;

	int       m_iMon;
	std::vector<Vec> m_solVector;

	bool m_bStoreVec;

	// for initial accn solve
	Mat       m_matMass;
	KSP       m_AccnKSP;

	// Matrix free ?
	bool  m_bMatrixFree;
};


int newmark::init() {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	int matsize;

	m_dBeta = 0.25; m_dGamma = 0.5;

	// Allocate memory for working vectors
	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vecSolution));
	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vecVelocity ));
	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vecAccn ));

	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vecRHS));

	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vec_du));
	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vec_dv));
	CHKERRQ(VecDuplicate(m_vecInitialSolution, &m_vec_da));

	CHKERRQ(VecGetLocalSize(m_vecInitialSolution, &matsize));

	if (m_bMatrixFree) {
		CHKERRQ(MatCreateShell(PETSC_COMM_WORLD, matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE, this, &m_matJacobian));
		CHKERRQ(MatShellSetOperation(m_matJacobian, MATOP_MULT, (void(*)(void))(MatMult)));
		CHKERRQ(MatShellSetOperation(m_matJacobian, MATOP_GET_DIAGONAL, (void(*)(void))(MatGetDiagonal)));
		// MatShell for initial accn solve ...
		CHKERRQ(MatCreateShell(PETSC_COMM_WORLD, matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE, this, &m_matMass));
		CHKERRQ(MatShellSetOperation(m_matMass, MATOP_MULT, (void(*)(void))(InitMatMult)));
		CHKERRQ(MatShellSetOperation(m_matMass, MATOP_GET_DIAGONAL, (void(*)(void))(InitMatGetDiagonal)));
	} else {
		// std::cout << "Assembling Matrices" << std::endl;
		double dt = m_ti->step;
		Mat M, K;
		// Get the matrices ...
		m_Mass->GetAssembledMatrix(&M, MATAIJ);
		// std::cout << "Done Mass" << std::endl;
		m_Stiffness->GetAssembledMatrix(&K, MATAIJ);

		// std::cout << "Done assembling Matrices" << std::endl;

		// Create the actual ones ...
		MatDuplicate(M, MAT_COPY_VALUES, &m_matMass);

		MatScale(M, 1.0/(m_dBeta*dt*dt));
		MatDuplicate(M, MAT_COPY_VALUES, &m_matJacobian);
		MatAXPY(m_matJacobian, -1.0, K, SAME_NONZERO_PATTERN);

		int npes;
		MPI_Comm_size(MPI_COMM_WORLD, &npes);

		if (npes == 1) {
			PetscTruth isSym;
			MatIsSymmetric(m_matJacobian, 1e-16, &isSym);
			if (isSym == PETSC_TRUE)
				std::cout << "Matrix is symmetric" << std::endl;
			else
				std::cout	<< "Matrix is not symmetric" <<	std::endl;
			// Write out the matrix ...
			/*
					PetscViewer        viewer;
					PetscViewerASCIIOpen(MPI_COMM_WORLD, "mat.m" , &viewer);
	
					PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
					MatView(m_matJacobian, viewer);
					// MatView(M, viewer);
					PetscViewerFlush(viewer);
					PetscViewerDestroy(viewer);
				*/
		}

		// Destroy the rest ...
		MatDestroy(K);
		MatDestroy(M);

		// Write out the matrix to a matlab file. 
		/*
		PetscViewer viewer;
		char viewername[100];
		sprintf(viewername,"lvMat.m");
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,viewername,&viewer);
		PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
		MatView(m_matJacobian, viewer);
		PetscViewerDestroy(viewer);
		std::cout << "Finished writing matrix" << std::endl;
		*/
	}

	// Create a KSP context to solve  @ every timestep
	CHKERRQ(KSPCreate(PETSC_COMM_WORLD, &m_ksp));
	CHKERRQ(KSPSetOptionsPrefix(m_ksp,"fwd_"));
	CHKERRQ(KSPSetOperators(m_ksp, m_matJacobian, m_matJacobian, SAME_NONZERO_PATTERN));
	CHKERRQ(KSPSetType(m_ksp,KSPCG));

	CHKERRQ(KSPSetFromOptions(m_ksp));

	// KSP for initial accn solve ...

	// Create a KSP context to solve  @ every timestep
	CHKERRQ(KSPCreate(PETSC_COMM_WORLD, &m_AccnKSP));
	CHKERRQ(KSPSetOptionsPrefix(m_AccnKSP,"mass_"));
	CHKERRQ(KSPSetOperators(m_AccnKSP, m_matMass, m_matMass, SAME_NONZERO_PATTERN));
	CHKERRQ(KSPSetType(m_AccnKSP,KSPCG));

	CHKERRQ(KSPSetFromOptions(m_AccnKSP));



	// initialize accn ...
	// a  = inv(M) * (P(:,1) - C * v(:,1) - K * u(:,1));
	// for now this is zero ...
	// CHKERRQ( VecZeroEntries( m_vecAccn ) );
	// Will be done before the solve ...
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return(0);
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_Solve"
int newmark::solve() {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	// Time info
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	m_dBeta = 0.25; m_dGamma = 0.5;
	// std::cout << "Newmark beta = " << m_dBeta << " and Gamma is " << m_dGamma << std::endl;

	unsigned NT = (int)(ceil((m_ti->stop - m_ti->start)/m_ti->step));

	// std::cout << RED"Number of Timesteps is "NRM << NT << std::endl;
	double temprtol;
	CHKERRQ(KSPGetTolerances(m_ksp, &temprtol,0,0,0));

	// Set initial conditions
	CHKERRQ( VecCopy( m_vecInitialSolution, m_vecSolution ) );
	CHKERRQ( VecCopy( m_vecInitialVelocity, m_vecVelocity ) );

	// Calculate the initial acceleration ...
	CHKERRQ( VecZeroEntries( m_vecAccn ) );
	//std::cout << "Setting RHS for initial accn" << std::endl;
	setAccnRHS();
	//std::cout << "Solving for initial accn" << std::endl;
	CHKERRQ( KSPSolve (m_AccnKSP, m_vecRHS, m_vecAccn) );
	// std::cout << "Done solving for initial accn" << std::endl;


	// Time stepping
	if ( !m_bIsAdjoint) {	// FORWARD PROBLEM
		// std::cout << RED"STARTING FORWARD SOLVE"NRM << std::endl;
		m_ti->current = m_ti->start;
		m_ti->currentstep = 0;
		monitor();
		while (m_ti->currentstep < NT) {
			// std::cout << "Step is " << m_ti->currentstep << std::endl;
			// std::cout << GRN"At time "RED << m_ti->current << NRM" ";
			double stime = MPI_Wtime();

			m_ti->currentstep++;
			m_ti->current += m_ti->step;

#ifdef __DEBUG__
			if (!rank) {
				std::cout << GRN"Forward Time is "YLW << std::setw(6) << std::setprecision(3) << m_ti->current << NRM"\r" << std::flush;
			}
#endif

			// Get the Right hand side of the ksp solve using the current solution
			// std::cout << "Setting RHS" << std::endl;
			setRHS();
			// std::cout << "Finished setting RHS" << std::endl;

#ifdef __DEBUG__
			double norm;
			CHKERRQ(VecNorm(m_vecRHS, NORM_INFINITY,&norm));
			std::cout << GRN"norm of the RHS @ "NRM << m_ti->current  << " = " << norm << std::endl;
#endif

			// clear du, dv, and da
			CHKERRQ ( VecZeroEntries(m_vec_du) );
			CHKERRQ ( VecZeroEntries(m_vec_dv) );
			CHKERRQ ( VecZeroEntries(m_vec_da) );

			// Solve the ksp using the current rhs and non-zero initial guess
			// CHKERRQ( KSPSetInitialGuessNonzero( m_ksp, PETSC_TRUE) );
			// std::cout << "Starting KSP Solve" << std::endl;
			CHKERRQ( KSPSolve (m_ksp, m_vecRHS, m_vec_du) );

			// std::cout << "Finished KSP Solve" << std::endl;

#ifdef __DEBUG__      
			double norm3;
			CHKERRQ ( VecNorm(m_vec_du, NORM_INFINITY, &norm3) );
			std::cout << "norm of the DU @ " << m_ti->current  << " = " << norm3 << std::endl;
#endif      

			double dt = m_ti->step;

			// compute dv and da
			CHKERRQ ( VecAXPY(m_vec_dv, m_dGamma/(m_dBeta*dt), m_vec_du));
			CHKERRQ ( VecAXPY(m_vec_dv, -m_dGamma/m_dBeta, m_vecVelocity));
			CHKERRQ ( VecAXPY(m_vec_dv, dt*(1.0 - m_dGamma/(2.0*m_dBeta)), m_vecAccn));

#ifdef __DEBUG__      
			double norm4;
			CHKERRQ ( VecNorm(m_vec_dv, NORM_INFINITY, &norm4) );
			std::cout << "norm of the DV @ " << m_ti->current  << " = " << norm4 << std::endl;
#endif      


			CHKERRQ ( VecAXPY(m_vec_da, 1.0/(m_dBeta*dt*dt), m_vec_du));
			CHKERRQ ( VecAXPY(m_vec_da, -1.0/(m_dBeta*dt), m_vecVelocity));
			CHKERRQ ( VecAXPY(m_vec_da, -1.0/(2.0*m_dBeta), m_vecAccn));

			// Now update solution, velocity and acceleration.
			CHKERRQ ( VecAXPY(m_vecSolution, 1.0, m_vec_du) );
			CHKERRQ ( VecAXPY(m_vecVelocity, 1.0, m_vec_dv) );
			CHKERRQ ( VecAXPY(m_vecAccn,     1.0, m_vec_da) );

#ifdef __DEBUG__
			double norm1;
			CHKERRQ ( VecNorm(m_vecSolution, NORM_INFINITY, &norm1) );
			std::cout << GRN"norm of the solution @ "NRM << m_ti->current  << " = "RED << norm1 << NRM << std::endl;

			double norm2;
			CHKERRQ ( VecNorm(m_vecAccn, NORM_INFINITY, &norm2) );
			std::cout << "norm of the accn @ " << m_ti->current  << " = " << norm2 << std::endl;
#endif
			monitor();

			// std::cout << "Finished time step in time " << MPI_Wtime() - stime << std::endl;
		}
#ifdef __DEBUG__
		if (!rank)
			std::cout << GRN"Finished Forward Solve"NRM << std::endl;
#endif
	} else {	// ADJOINT PROBLEM
		m_ti->current = m_ti->start;
		m_ti->currentstep = 0;
		monitor();
		while (m_ti->currentstep < NT) {
			m_ti->currentstep++;
			m_ti->current += m_ti->step;

#ifdef __DEBUG__
			if (!rank) {
				std::cout << GRN"Adjoint Time is "YLW << std::setw(6) << std::setprecision(3) << m_ti->stop - m_ti->current << NRM"\r" << std::flush;
				// std::cout << GRN"Time is "YLW << m_ti->current << NRM"\t:\r";
			}
#endif

			// Get the right hand side of the ksp solve for the adjoint problem
			setRHS();

			// clear du, dv, and da
			CHKERRQ(VecZeroEntries(m_vec_du));
			CHKERRQ(VecZeroEntries(m_vec_dv));
			CHKERRQ(VecZeroEntries(m_vec_da));

			// Solve ksp using the current rhs and non-zero initial guess
			// CHKERRQ( KSPSetInitialGuessNonzero(m_ksp,PETSC_TRUE) );
			CHKERRQ( KSPSolve(m_ksp, m_vecRHS, m_vec_du) );

			double dt = m_ti->step;
			// compute dv and da
			CHKERRQ ( VecAXPY(m_vec_dv, m_dGamma/(m_dBeta*dt), m_vec_du));
			CHKERRQ ( VecAXPY(m_vec_dv, -m_dGamma/m_dBeta, m_vecVelocity));
			CHKERRQ ( VecAXPY(m_vec_dv, dt*(1.0 - m_dGamma/(2.0*m_dBeta)), m_vecAccn));

			CHKERRQ ( VecAXPY(m_vec_da, 1.0/(m_dBeta*dt*dt), m_vec_du));
			CHKERRQ ( VecAXPY(m_vec_da, -1.0/(m_dBeta*dt), m_vecVelocity));
			CHKERRQ ( VecAXPY(m_vec_da, -1.0/(2.0*m_dBeta), m_vecAccn));

			// Now update solution, velocity and acceleration.
			CHKERRQ ( VecAXPY(m_vecSolution, 1.0, m_vec_du) );
			CHKERRQ ( VecAXPY(m_vecVelocity, 1.0, m_vec_dv) );
			CHKERRQ ( VecAXPY(m_vecAccn,     1.0, m_vec_da) );

			monitor();
		}
#ifdef __DEBUG__
		if (!rank)
			std::cout << GRN"Finished Adjoint Solve"NRM << std::endl;
#endif
	}

	// std::cout << RED"Finished Solve"NRM << std::endl;
#ifdef __DEBUG__
	std::cout << "Leaving " << __func__ << std::endl;
#endif
	return(0);
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_MassMult"
void newmark::massMatMult(Vec In, Vec Out) {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	VecZeroEntries(Out); /* Clear to zeros*/
	//  std::cout << "dt is " << dt << std::endl;
	m_Mass->MatVec(In, Out, 1.0);
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_MassGetDiagonal"
void newmark::massGetDiagonal(Vec diag) {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	VecZeroEntries(diag) ; /* Clear to zeros*/
	m_Mass->MatGetDiagonal(diag, 1.0);
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_JacMatMult"
void newmark::jacobianMatMult(Vec In, Vec Out) {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	VecZeroEntries(Out); /* Clear to zeros*/
	double dt = m_ti->step;
	//  std::cout << "dt is " << dt << std::endl;
	m_Stiffness->MatVec(In, Out, -1.0);
	m_Mass->MatVec(In, Out, 1.0/(m_dBeta*dt*dt));
	// double norm1 = 0;
	// ierr = VecNorm(In, NORM_INFINITY,&norm1);
	// printf("norm of the Matvec after MassTerm %d\n", norm1);
	if (m_bDamp) {
		// std::cout << RED"Adding Damping"NRM << std::endl;
		m_Damping->MatVec(In, Out, m_dGamma/(m_dBeta*dt));
	}
	// ierr = VecNorm(In, NORM_INFINITY,&norm1);
	// std::cout << "norm of the Matvec after StiffnessTerm " << norm << std::endl;
	// printf("norm of the Matvec after StiffnessTerm %d\n", norm1);
	// norm1 = 0.0;
	// return;
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_JacMatGetDiag"
void newmark::jacobianGetDiagonal(Vec diag) {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	//VecSet(diag, 0.0); /* Clear to zeros*/
	VecZeroEntries(diag);

	double dt = m_ti->step;

	// double norm1;
	// PetscInt ierr;
	m_Mass->MatGetDiagonal(diag, 1.0/(m_dBeta*dt*dt));
	//m_Mass->MatGetDiagonal(diag, 1.0);
	// ierr = VecNorm(diag, NORM_2,&norm1);
	// printf("norm of the Diag after MassTerm %f\n", norm1);
	if (m_bDamp) {
		// std::cout << RED"Adding Damping"NRM << std::endl;
		m_Damping->MatGetDiagonal(diag, m_dGamma/(m_dBeta*dt));
	}
	m_Stiffness->MatGetDiagonal(diag, -1.0);
	// ierr = VecNorm(diag, NORM_2,&norm1);
	// printf("norm of the Diag after StiffTerm %f\n", norm1);

	// PetscScalar sum;
	// ierr = VecSum(diag, &sum);
	// printf("Trace of the Diag after StiffTerm %f\n", sum);


#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_SetAccnRHS"
bool newmark::setAccnRHS() {
	// std::cout << "Entering " << __func__ << std::endl;
	VecZeroEntries(m_vecRHS);

	unsigned NT = (int)(ceil((m_ti->stop - m_ti->start)/m_ti->step));
	if ( !m_bIsAdjoint) {	// FORWARD PROBLEM
		m_Force->addVec(m_vecRHS, 1.0, 0);	// Add force to right hand side
	} else {
		m_Force->addVec(m_vecRHS, 1.0, NT);	 // Add force to right hand side 
	}

	if (m_bDamp) {
		m_Damping->MatVec(m_vecInitialVelocity, m_vecRHS, -1.0 );
	}
	m_Stiffness->MatVec(m_vecInitialSolution, m_vecRHS, -1.0 );
	return true;
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_SetRHS"
bool newmark::setRHS() {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << " for timestep " << m_ti->currentstep << std::endl;
#endif  

	// double norm;  
	unsigned NT = (int)(ceil((m_ti->stop - m_ti->start)/m_ti->step));

	VecZeroEntries(m_vecRHS);

	m_Mass->MatVec(m_vecAccn, m_vecRHS, 1.0/(2.0*m_dBeta));

	if (m_bDamp) {
		m_Damping->MatVec(m_vecAccn, m_vecRHS, m_ti->step*(m_dGamma/(2*m_dBeta) -1) );
	}
	m_Mass->MatVec(m_vecVelocity, m_vecRHS, 1.0/(m_dBeta*m_ti->step));
	if (m_bDamp) {
		m_Damping->MatVec(m_vecVelocity, m_vecRHS, m_dGamma/m_dBeta );
	}

	if ( !m_bIsAdjoint) {	// FORWARD PROBLEM
		// if (m_ti->currentstep == 1) {
		m_Force->addVec(m_vecRHS, 1.0, m_ti->currentstep);	/* Add force to right hand side*/
		m_Force->addVec(m_vecRHS, -1.0, m_ti->currentstep - 1);	 /* Add force to right hand side*/
		// } else {
		//   m_Force->addVec(m_vecRHS, 0.5, m_ti->currentstep);  /* Add force to right hand side*/
		//   m_Force->addVec(m_vecRHS, -0.5, m_ti->currentstep - 2);  /* Add force to right hand side*/
		// }
	} else {
		// if (m_ti->currentstep == 1) {
		m_Force->addVec(m_vecRHS, 1.0, NT - m_ti->currentstep);	 /* Add force to right hand side*/
		m_Force->addVec(m_vecRHS, -1.0, NT - m_ti->currentstep + 1);	/* Add force to right hand side*/
		// } else {
		//  m_Force->addVec(m_vecRHS, 0.5, NT - m_ti->currentstep);  /* Add force to right hand side*/
		//   m_Force->addVec(m_vecRHS, -0.5, NT - m_ti->currentstep + 2);  /* Add force to right hand side*/
		// }
	}

#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

#undef __FUNCT__
#define __FUNCT__ "Newmark_Monitor"
int newmark::monitor() {
#ifdef __DEBUG__  
	std::cout << RED"In monitor"NRM << std::endl;
	std::cout << RED"imon is "NRM << m_iMon << std::endl; 
#endif

	if (fmod(m_ti->currentstep,(double)(m_iMon)) < 0.0001) {
		// double norm;
		int ierr;
		Vec tempSol;
		ierr = VecDuplicate(m_vecSolution,&tempSol);  CHKERRQ(ierr);
#ifdef __DEBUG__
		//		VecNorm(m_vecSolution,NORM_INFINITY,&norm);
		//		PetscPrintf(0,"solution norm b4 push back %f\n",norm);
#endif

		// m_bStoreVec=false;


		if (m_bStoreVec) {
			ierr = VecCopy(m_vecSolution,tempSol); CHKERRQ(ierr);
			// std::cout << YLW"Pushing to solution Vector"NRM << std::endl;
			if ( !m_bIsAdjoint) {	// FORWARD PROBLEM
				m_solVector.push_back(tempSol);
			} else {
				m_solVector.insert(m_solVector.begin(), tempSol);
			}
		} else {
			// will write out the file as,
			// Def.nt.rank.raw



			PetscScalar * buff;
			int sz;
			VecGetSize(m_vecSolution, &sz);
			CHKERRQ ( VecGetArray(m_vecSolution, &buff) );

			// int rank;
			// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			/*
			double norm;
			VecNorm(m_vecSolution, NORM_INFINITY, &norm);
			std::cout << "Norm at step " << m_ti->currentstep << " is " << norm << std::endl;
			*/
			/*
			float *fbuff = new float[sz];
			for (int i=0; i<sz; i++) {
				fbuff[i] = (float)buff[i];
			}
			*/
			char fname[256];
			char hdrname[256];
			sprintf(fname, "Def.%.3d.raw", m_ti->currentstep);
			// sprintf(fname, "Data_%.4d.%.4d.raw", m_ti->currentstep, rank);
			std::ofstream out(fname, std::ios::binary);
			out.write((char *)buff, sizeof(double)*sz );
			out.close();

			// MHD header
			int Ns = (int) round(pow(sz/3, 0.33333));
			sprintf(hdrname, "Def.%.3d.mhd", m_ti->currentstep);
			out.open(hdrname);
			out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
			out << "ElementSpacing =  2.73067 2.73067 2.73067" << std::endl;
			out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
			out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_DOUBLE" << std::endl;
			out << "ElementDataFile = " << fname << std::endl;
			out.close();

			// delete [] fbuff;
			CHKERRQ ( VecRestoreArray(m_vecSolution, &buff) );

			if (m_ti->currentstep > 1) {
				// write out the force ...
				/*
				double norm;
				VecNorm(m_vecRHS, NORM_INFINITY, &norm);
				std::cout << "Force Norm at step " << m_ti->currentstep << " is " << norm << std::endl;
				*/

				CHKERRQ ( VecGetArray(m_vecRHS, &buff) );
				sprintf(fname, "force.%.3d.raw", m_ti->currentstep);
				out.open(fname, std::ios::binary);
				out.write((char *)buff, sizeof(double)*sz );
				out.close();

				// MHD header
				int Ns = (int)round(pow(sz/3, 0.33333));
				sprintf(hdrname, "force.%.3d.mhd", m_ti->currentstep);
				out.open(hdrname);
				out << "ObjectType = Image" << std::endl << "NDims = 3" << std::endl << "BinaryData = True" << std::endl << "BinaryDataByteOrderMSB = False" << std::endl << "Offset = 0 0 0" << std::endl;
				out << "ElementSpacing =  2.73067 2.73067 2.73067" << std::endl;
				out << "DimSize = " << Ns << " " << Ns << " " << Ns << std::endl;
				out << "ElementNumberOfChannels = 3" << std::endl << "ElementType = MET_DOUBLE" << std::endl;
				out << "ElementDataFile = " << fname << std::endl;
				out.close();

				// delete [] fbuff;
				CHKERRQ ( VecRestoreArray(m_vecRHS, &buff) );
			}
		}

#ifdef __DEBUG__
		VecNorm(m_solVector[m_solVector.size()-1],NORM_INFINITY,&norm);
		PetscPrintf(0,"solution norm after push back %f\n",norm);
#endif
	}
	return(0);
}

#endif

