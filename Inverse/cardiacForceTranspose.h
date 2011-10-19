/**
 *  @file	cardiacForceTranspose.h
 *  @brief	Main class to get the Transpose of the activation
 *         force within the myocardium. This is required for the
 *         gradient computation.
 *  @author Hari Sundar
 *  @date	  5/15/7
 * 
 *  Main class to get the Transpose of the activation force
 *         within the myocardium. This is required for the
 *         gradient computation.
 */

#ifndef __CARDIAC_FORCE_TRANSPOSE_H_
#define __CARDIAC_FORCE_TRANSPOSE_H_

#include <vector>

#include "feVector.h"

/**
 *  @brief	Main class to get the Transpose of the activation
 *         force within the myocardium. This is required for the
 *         gradient computation.
 *  @author Hari Sundar
 *  @date	  6/7/7
 * 
 *  Main class to get the Transpose of the activation
 *         force within the myocardium. This is required for the
 *         gradient computation. The fiber orientations within the myocardium
 *  need to be specified. In addition the scalar
 *  activation(\f$\tau\f$) term needs to be specified a priori.
 */
class cardiacForceTranspose : public feVector<cardiacForceTranspose> {
public: 
	cardiacForceTranspose(daType da);

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
	void setAdjoints(std::vector<Vec> adjoints) {
		m_vecAdjoints = adjoints;
	}

	void setFiberOrientations(Vec fib) {
		fibersVec = fib;
	}

	void setDA3d(DA da) {
		m_da3D = da;
	}

private:
	std::vector<Vec>         m_vecAdjoints;
	Vec                      fibersVec;

	void *                    m_lambda;
	void *                    m_fibers;

	double     m_dHx;

	DA m_da3D;

	double xFac, yFac, zFac;
	unsigned int maxD;
};

cardiacForceTranspose::cardiacForceTranspose(daType da) {
#ifdef __DEBUG__
	assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
	m_daType = da;
	m_DA    = NULL;
	m_octDA   = NULL;
	m_stencil = NULL;

	m_lambda = NULL;
	m_fibers = NULL;

	// initialize the stencils ...
	initStencils();
	if (da == OCT)
		initOctLut();

}
bool cardiacForceTranspose::initStencils() {
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

bool cardiacForceTranspose::preAddVec() {
#ifdef __DEBUG__    
	std::cout << "Entering cFT::" << __func__ << std::endl;
#endif    
	if (m_daType == PETSC) {
		// compute Hx
		PetscInt mx,my,mz;
		PetscScalar ***lam;	// 
		PetscScalar ***fibers; // 

		// std::cout << "Getting adjoints" << std::endl;
		CHKERRQ( DAVecGetArray(m_da3D, m_vecAdjoints[m_iCurrentDynamicIndex], &lam) );
		// std::cout << "Getting fibers" << std::endl;
		CHKERRQ( DAVecGetArray(m_da3D, fibersVec, &fibers) );
		// std::cout << "Got it " << std::endl;
		m_lambda = lam;
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
		m_octDA->vecGetBuffer(m_vecAdjoints[m_iCurrentDynamicIndex], tau, false, false, true, m_uiDof);
		m_octDA->vecGetBuffer(fibersVec, fibers, false, false, true, m_uiDof);
		m_lambda = tau;
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

bool cardiacForceTranspose::postAddVec() {
#ifdef __DEBUG__
	std::cout << "Entering " << __func__ << std::endl;
#endif
	if ( m_daType == PETSC) {
		PetscScalar ***tau = (PetscScalar ***)m_lambda;
		CHKERRQ( DAVecRestoreArray(m_da3D, m_vecAdjoints[m_iCurrentDynamicIndex], &tau) );

		PetscScalar ***fibers = (PetscScalar ***)m_fibers;
		CHKERRQ( DAVecRestoreArray(m_da3D, fibersVec, &fibers) );
	} else {
		PetscScalar *tau = (PetscScalar *)m_lambda;
		m_octDA->vecRestoreBuffer(m_vecAdjoints[m_iCurrentDynamicIndex], tau, false, false, true, m_uiDof);
		PetscScalar *fibers = (PetscScalar *)m_fibers;
		m_octDA->vecRestoreBuffer(fibersVec, fibers, false, false, true, m_uiDof);
	}
#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

bool cardiacForceTranspose::ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale) {
#ifdef __DEBUG__
	std::cout << "Entering " << __func__ << " " << i << ", " << j << ", " << k << std::endl;
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

	PetscScalar ***lambda     = (PetscScalar ***) m_lambda;
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
		nx = ny = nz = 0.0;
	}

	double gix=0.0, giy=0.0, giz=0.0;
	// first compute q = (n*n') p
	double * qi = new double [8*dof];

	for (int q = 0; q < 8; q++) {
		gix = lambda[idx[q][0]][idx[q][1]][dof*idx[q][2]];
		giy = lambda[idx[q][0]][idx[q][1]][dof*idx[q][2]+1];
		giz = lambda[idx[q][0]][idx[q][1]][dof*idx[q][2]+1];

		qi[3*q]   = gix*nx*nx + giy*nx*ny + giz*nx*nz;
		qi[3*q+1] = gix*nx*ny + giy*ny*ny + giz*ny*nz;
		qi[3*q+2] = gix*nx*nz + giy*nz*ny + giz*nz*nz;
	}



	for (int q = 0; q < 8; q++) {
		// now compute g = A' q;
		// which is 
		// g(i)  = \sum_j A(j,i) \tau(j)
		gix = 0.0;
		for (int r = 0; r < 8; r++) {
			gix += Ajk[3*r][q] * qi[3*r]  +  
						 Ajk[3*r+1][q] * qi[3*r+1] +
						 Ajk[3*r+2][q] * qi[3*r+2];
		}        
		in[idx[q][0]][idx[q][1]][idx[q][2]] += gix*stencilScale;      
	}            

	delete [] qi;
#ifdef __DEBUG__    
	std::cout << "Leaving " << __func__ << std::endl;
#endif    
	return true;
}

bool cardiacForceTranspose::ElementalAddVec(unsigned int index, PetscScalar *in, double scale) {
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

	PetscScalar *lambda     = (PetscScalar *) m_lambda;
	PetscScalar *fibers  = (PetscScalar *) m_fibers;

	double nx,ny,nz; //, nx0, ny0, nz0, snx, sny, snz;
	nx = fibers[m_uiDof*index];
	ny = fibers[m_uiDof*index + 1];
	nz = fibers[m_uiDof*index + 2];
	if ( sqrt(nx*nx+ ny*ny + nz*nz) < 0.001) {
		nx = ny = nz = 0.0;
	}

	double gix=0.0, giy=0.0, giz=0.0;
	// first compute q = (n*n') p
	double * qi = new double [8*m_uiDof];

	for (int q = 0; q < 8; q++) {
		gix = lambda[m_uiDof * idx[q]];
		giy = lambda[m_uiDof * idx[q] + 1];
		giz = lambda[m_uiDof * idx[q] + 2];

		qi[3*q]   = gix*nx*nx + giy*nx*ny + giz*nx*nz;
		qi[3*q+1] = gix*nx*ny + giy*ny*ny + giz*ny*nz;
		qi[3*q+2] = gix*nx*nz + giy*nz*ny + giz*nz*nz;
	}

	for (int q = 0; q < 8; q++) {
		// now compute g = A' q;
		// which is 
		// g(i)  = \sum_j A(j,i) \tau(j)
		gix = 0.0;
		for (int r = 0; r < 8; r++) {
			gix += A[8*3*r+q] * qi[3*r]  +  
						 A[8*3*(r+1)+q] * qi[3*r+1] +
						 A[8*3*(r+2)+q] * qi[3*r+2];
		}        
		in[idx[q]] += gix*stencilScale;      
	}          

#ifdef __DEBUG__  
	std::cout << "Leaving " << __func__ << std::endl;
#endif  
	return true;
}

#endif
