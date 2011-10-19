/**
 *  @file	cardiacForce.h
 *  @brief	Main class to set the activation force within the
 *         myocardium.
 *  @author Hari Sundar
 *  @date	  5/15/7
 * 
 *  Main  class to set  the activation force within the
 *  myocardium.
 */

#ifndef __CARDIAC_FORCE_H_
#define __CARDIAC_FORCE_H_

#include <vector>

#include "feVector.h"

/**
 *  @brief	Main class to set the activation force within the
 *         myocardium.
 *  @author Hari Sundar
 *  @date	  6/7/7
 * 
 *  Main  class to set  the activation force within the
 *  myocardium. The fiber orientations within the myocardium
 *  need to be specified. In addition the scalar
 *  activation(\f$\tau\f$) term needs to be specified a priori.
 */
class cardiacForce : public feVector<cardiacForce> {
public: 
	cardiacForce(daType da);

	bool initStencils();

	inline bool ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale);
	inline bool ElementalAddVec(unsigned int index, PetscScalar *in, double scale);

	bool preAddVec();
	bool postAddVec();

	/**
	 * Set the activation forces.
	 * 
	 * @param actVec
	 */
	void setActivationVec(std::vector<Vec> actVec) {
		tauVec = actVec;
	}

	void setFiberOrientations(Vec fib) {
		fibersVec = fib;
	}

	Vec getFiberOrientations() {
		return fibersVec;
	}

private:
	std::vector<Vec>         tauVec;
	Vec                      fibersVec;

	void *                    m_tau;
	void *                    m_fibers;

	double     m_dHx;

	double xFac, yFac, zFac;
	unsigned int maxD;
};

cardiacForce::cardiacForce(daType da) {
#ifdef __DEBUG__
	assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
	m_daType = da;
	m_DA    = NULL;
	m_octDA   = NULL;
	m_stencil = NULL;

	m_tau = NULL;
	m_fibers = NULL;

	// initialize the stencils ...
	initStencils();
	if (da == OCT)
		initOctLut();

}
bool cardiacForce::initStencils() {
#ifdef __DEBUG__    
	std::cout << "Entering " << __func__ << std::endl;
#endif    
	typedef double* doublePtr;

	if (m_daType == PETSC) {
		double Bjk[24][8] =   {
			{-0.222222222222222,-0.222222222222222,-0.111111111111111,-0.111111111111111,-0.111111111111111,-0.111111111111111,-0.055555555555556,-0.055555555555556},
			{-0.222222222222222,-0.111111111111111,-0.222222222222222,-0.111111111111111,-0.111111111111111,-0.055555555555556,-0.111111111111111,-0.055555555555556},
			{-0.222222222222222,-0.111111111111111,-0.111111111111111,-0.055555555555556,-0.222222222222222,-0.111111111111111,-0.111111111111111,-0.055555555555556},
			{0.222222222222222,0.222222222222222,0.111111111111111,0.111111111111111,0.111111111111111,0.111111111111111,0.055555555555556,0.055555555555556},
			{-0.111111111111111,-0.222222222222222,-0.111111111111111,-0.222222222222222,-0.055555555555556,-0.111111111111111,-0.055555555555556,-0.111111111111111},
			{-0.111111111111111,-0.222222222222222,-0.055555555555556,-0.111111111111111,-0.111111111111111,-0.222222222222222,-0.055555555555556,-0.111111111111111},
			{-0.111111111111111,-0.111111111111111,-0.222222222222222,-0.222222222222222,-0.055555555555556,-0.055555555555556,-0.111111111111111,-0.111111111111111},
			{0.222222222222222,0.111111111111111,0.222222222222222,0.111111111111111,0.111111111111111,0.055555555555556,0.111111111111111,0.055555555555556},
			{-0.111111111111111,-0.055555555555556,-0.222222222222222,-0.111111111111111,-0.111111111111111,-0.055555555555556,-0.222222222222222,-0.111111111111111},
			{0.111111111111111,0.111111111111111,0.222222222222222,0.222222222222222,0.055555555555556,0.055555555555556,0.111111111111111,0.111111111111111},
			{0.111111111111111,0.222222222222222,0.111111111111111,0.222222222222222,0.055555555555556,0.111111111111111,0.055555555555556,0.111111111111111},
			{-0.055555555555556,-0.111111111111111,-0.111111111111111,-0.222222222222222,-0.055555555555556,-0.111111111111111,-0.111111111111111,-0.222222222222222},
			{-0.111111111111111,-0.111111111111111,-0.055555555555556,-0.055555555555556,-0.222222222222222,-0.222222222222222,-0.111111111111111,-0.111111111111111},
			{-0.111111111111111,-0.055555555555556,-0.111111111111111,-0.055555555555556,-0.222222222222222,-0.111111111111111,-0.222222222222222,-0.111111111111111},
			{0.222222222222222,0.111111111111111,0.111111111111111,0.055555555555556,0.222222222222222,0.111111111111111,0.111111111111111,0.055555555555556},
			{0.111111111111111,0.111111111111111,0.055555555555556,0.055555555555556,0.222222222222222,0.222222222222222,0.111111111111111,0.111111111111111},
			{-0.055555555555556,-0.111111111111111,-0.055555555555556,-0.111111111111111,-0.111111111111111,-0.222222222222222,-0.111111111111111,-0.222222222222222},
			{0.111111111111111,0.222222222222222,0.055555555555556,0.111111111111111,0.111111111111111,0.222222222222222,0.055555555555556,0.111111111111111},
			{-0.055555555555556,-0.055555555555556,-0.111111111111111,-0.111111111111111,-0.111111111111111,-0.111111111111111,-0.222222222222222,-0.222222222222222},
			{0.111111111111111,0.055555555555556,0.111111111111111,0.055555555555556,0.222222222222222,0.111111111111111,0.222222222222222,0.111111111111111},
			{0.111111111111111,0.055555555555556,0.222222222222222,0.111111111111111,0.111111111111111,0.055555555555556,0.222222222222222,0.111111111111111},
			{0.055555555555556,0.055555555555556,0.111111111111111,0.111111111111111,0.111111111111111,0.111111111111111,0.222222222222222,0.222222222222222},
			{0.055555555555556,0.111111111111111,0.055555555555556,0.111111111111111,0.111111111111111,0.222222222222222,0.111111111111111,0.222222222222222},
			{0.055555555555556,0.111111111111111,0.111111111111111,0.222222222222222,0.055555555555556,0.111111111111111,0.111111111111111,0.222222222222222}
		};

		double** Ajk = new doublePtr[24];
		for (int j=0;j<24;j++) {
			Ajk[j] = new double[8];
			for (int k=0;k<8;k++) {
				Ajk[j][k] = Bjk[j][k];
			}//end k
		}//end j
		m_stencil = Ajk;
	} else {
		m_stencil = NULL;
		// allocate memory for the stencils 
		double *K = new double[8*18*24*8];
		// read in the stencils from the file ...
		std::ifstream in("FiberForce.inp");
		in.read((char *)K, 8*18*24*8*8);
		in.close();
		m_stencil = K;
	}
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

bool cardiacForce::preAddVec() {
#ifdef __DEBUG__    
	std::cout << "Entering " << __func__ << std::endl;
#endif    
	if (m_daType == PETSC) {
		// compute Hx
		PetscInt mx,my,mz;
		PetscScalar *tau;	// 
		PetscScalar ***fibers; // 

		double tauNorm;
		VecNorm(tauVec[m_iCurrentDynamicIndex], NORM_2, &tauNorm);
		// tauNorm = tauNorm/pow(Ns,1.5);
		// std::cout << "Activation Norm is " << tauNorm << std::endl;


		CHKERRQ( VecGetArray(tauVec[m_iCurrentDynamicIndex], &tau) );
		CHKERRQ( DAVecGetArray(m_DA, fibersVec, &fibers) );
		m_tau = tau;
		m_fibers = fibers;

		CHKERRQ( DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ); 

		m_dHx = m_dLx/(mx -1);
		m_dHx = m_dHx*m_dHx; //*m_dHx;
		m_dHx /= 4.0;
	} else {
		unsigned int maxD = m_octDA->getMaxDepth();
		double xFac,yFac,zFac;

		PetscScalar *tau; 
		PetscScalar *fibers; 
		// Get arrays
		m_octDA->vecGetBuffer(tauVec[m_iCurrentDynamicIndex], tau, false, false, true, 1);
		m_octDA->vecGetBuffer(fibersVec, fibers, false, false, true, m_uiDof);
		m_tau = tau;
		m_fibers = fibers;

		// Get the  x,y,z factors 
		xFac = 1.0/((double)(1<<(maxD-1)));
		if (m_octDA->getDimension() > 1) {
			yFac = 1.0/((double)(1<<(maxD-1)));
			if (m_octDA->getDimension() > 2) {
				zFac = 1.0/((double)(1<<(maxD-1)));
			}
		}
	}
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

bool cardiacForce::postAddVec() {
#ifdef __DEBUG__
	std::cout << "Entering " << __func__ << std::endl;
#endif
	if ( m_daType == PETSC) {
		PetscScalar *tau = (PetscScalar *)m_tau;
		CHKERRQ( VecRestoreArray(tauVec[m_iCurrentDynamicIndex], &tau) );

		PetscScalar ***fibers = (PetscScalar ***)m_fibers;
		CHKERRQ( DAVecRestoreArray(m_DA, fibersVec, &fibers) );
	} else {
		PetscScalar *tau = (PetscScalar *)m_tau;
		m_octDA->vecRestoreBuffer(tauVec[m_iCurrentDynamicIndex], tau, false, false, true, 1);
		PetscScalar *fibers = (PetscScalar *)m_fibers;
		m_octDA->vecRestoreBuffer(fibersVec, fibers, false, false, true, m_uiDof);
	}
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

#undef __FUNCT__
#define __FUNCT__ "cardiacForce_ElementalAddVec"
bool cardiacForce::ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale) {
#ifdef __DEBUG__
	std::cout << "Entering " << __func__ << std::endl;
#endif
	int dof=3;
	int idx[8][3]={
		{k, j, i},
		{k,j,i+1},
		{k,j+1,i},
		{k,j+1,i+1},
		{k+1, j, i},
		{k+1,j,i+1},
		{k+1,j+1,i},
		{k+1,j+1,i+1}               
	};             

	PetscScalar *tau     = (PetscScalar *) m_tau;
	PetscScalar ***fibers  = (PetscScalar ***) m_fibers;

	int x,y,z,m,n,p;
	CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

	double stencilScale =  m_dHx*scale;
	double **Ajk = (double **)m_stencil;

	// g = A\tau
	// f = (n.g)n
	double nx,ny,nz; //, nx0, ny0, nz0, snx, sny, snz;
	nx = fibers[k][j][dof*i];
	ny = fibers[k][j][dof*i + 1];
	nz = fibers[k][j][dof*i + 2];
	if ( sqrt(nx*nx+ ny*ny + nz*nz) < 0.001) {
		nx = ny = nz = 0.01;
	}
	/*
	nx0=nx;ny0=ny;nz0=nz;
	snx=nx;sny=ny;snz=nz;
	for (int q=1; q<8;q++) {
		nx = fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]];
		ny = fibers[idx[q][0]][idx[q][1]][dof*idx[q][2] + 1];
		nz = fibers[idx[q][0]][idx[q][1]][dof*idx[q][2] + 2];
		if ( (nx*nx0 + ny*ny0 + nz*nz0) < 0 ) {
			snx -= nx; sny -= ny; snz -= nz;  
		} else {
			snx += nx; sny += ny; snz += nz;  
		}
	}
	nx = snx/8; ny = sny/8; nz = snz/8;
	*/
	for (int q = 0; q < 8; q++) {
		double gix=0.0, giy=0.0, giz=0.0;
		double qix=0.0, qiy=0.0, qiz=0.0;

		// first compute g = A\tau
		// which is 
		// g(j)  = \sum_i A(j,i) \tau(i)
		for (int r = 0; r < 8; r++) {
			PetscScalar _tau = tau [ ((idx[r][0]-z)*n + idx[r][1] - y)*m + idx[r][2]-x ];	// tau[idx[r][0]][idx[r][1]][idx[r][2]];
			_tau = 10.0*_tau;

			gix += Ajk[3*q][r]   * _tau; 
			giy += Ajk[3*q+1][r] * _tau;
			giz += Ajk[3*q+2][r] * _tau;
		}        

		// fac =  <n,g>
		/*
		double fac = gix*fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]] + 
								 giy*fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+1] + 
								 giz*fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+2];						 
		*/

		// q = (n*n')g 
		qix = gix*nx*nx + giy*nx*ny + giz*nx*nz;
		qiy = gix*nx*ny + giy*ny*ny + giz*ny*nz;
		qiz = gix*nx*nz + giy*nz*ny + giz*nz*nz;
		// qix=1.0, qiy=1.0, qiz=1.0;

		in[idx[q][0]][idx[q][1]][dof*idx[q][2]]   +=  qix*stencilScale;	//* fac * fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]];  
		in[idx[q][0]][idx[q][1]][dof*idx[q][2]+1] +=  qiy*stencilScale;	//* fac * fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+1];
		in[idx[q][0]][idx[q][1]][dof*idx[q][2]+2] +=  qiz*stencilScale;	//* fac * fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+2];
	}                        
#ifdef __DEBUG__    
	std::cout << "Leaving " << __func__ << std::endl;
#endif    
	return true;
}

bool cardiacForce::ElementalAddVec(unsigned int index, PetscScalar *in, double scale) {
#ifdef __DEBUG__  
	std::cout << "Entering " << __func__ << std::endl;
#endif  
	unsigned int lev = m_octDA->getLevel(index);
	double hx = xFac*(1<<(maxD - lev));
	// double hy = yFac*(1<<(maxD - lev));
	// double hz = zFac*(1<<(maxD - lev));

	double stencilScale = scale*hx*hx/4.0;

	// stdElemType elemType;
	ot::DA::index idx[8];

	unsigned int chNum = m_octDA->getChildNumber();
	m_octDA->getNodeIndices(idx); 
	unsigned char hangingMask = m_octDA->getHangingNodeIndex( m_octDA->curr() );    
	unsigned char eType = getEtype(hangingMask, chNum);

	// get A which are the correct 24x8 matrices.
	double *K = (double *)m_stencil;
	double *A = K + (chNum*18 + eType)*24*8;

	PetscScalar *tau     = (PetscScalar *) m_tau;
	PetscScalar *fibers  = (PetscScalar *) m_fibers;

	// g = A\tau
	// f = (n.g)n

	for (int k = 0;k < 8; k++) {
		double gix=0.0, giy=0.0, giz=0.0;
		// first compute g = A\tau
		// which is 
		// g(k)  = \sum_j A(k,j) \tau(j)
		for (int j=0;j < 8; j++) {
			gix += A[8*3*k+j]   * tau[idx[j]];
			giy += A[8*(3*k+1)+j] * tau[idx[j]];
			giz += A[8*(3*k+2)+j] * tau[idx[j]];
		}//end for j

		// fac =  <n,g>
		double fac = gix*fibers[m_uiDof*idx[k]] + 
								 giy*fibers[m_uiDof*idx[k]+1] + 
								 giz*fibers[m_uiDof*idx[k]+2];

		in[m_uiDof*idx[k]]   += fac * stencilScale * fibers[m_uiDof*idx[k]];   
		in[m_uiDof*idx[k]+1] += fac * stencilScale * fibers[m_uiDof*idx[k]+1]; 
		in[m_uiDof*idx[k]+2] += fac * stencilScale * fibers[m_uiDof*idx[k]+2]; 
	}//end for k                      

#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

#endif
