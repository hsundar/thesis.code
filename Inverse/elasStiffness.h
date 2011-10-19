/**
 *  @file	elasStiffness.h
 *  @brief	Main class for finite element assembly of a linear elasticity
 *           stiffness matrix.
 *  @author	Hari Sundar
 *  @date	5/8/7
 *
 *  Main class for finite element assembly of linear elasticity stiffness matrix. The Lame 
 *  parameter images, i.e., the material properties, \f$\mu\f$ and \f$\lambda\f$, need to
 *  be specified. MatVec functions are virtual and in most cases need not be 
 *  specified as they are obtained from the parent feMatrix class. In most cases 
 *  only the assembly of the element matrix needs to be  done within this or a 
 *  derived class.
 */
#ifndef _ELAS_STIFFNESS_H_
#define _ELAS_STIFFNESS_H_

#include "feMatrix.h"
#include <fstream>

/**
 *  @brief	Main class for finite element assembly of a linear elasticity
 *           stiffness matrix.
 *  @author	Hari Sundar
 *  @date	5/8/7
 *
 *  Main class for finite element assembly of linear elasticity stiffness matrix. The Lame 
 *  parameter images, i.e., the material properties, \f$\mu\f$ and \f$\lambda\f$, need to
 *  be specified. MatVec functions are virtual and in most cases need not be 
 *  specified as they are obtained from the parent feMatrix class. In most cases 
 *  only the assembly of the element matrix needs to be  done within this or a 
 *  derived class.
 */

class elasStiffness : public feMatrix<elasStiffness> {
public:

  enum exhaustiveElemType {
    //Order: 654321
    //YZ ZX Z XY Y X
    ET_N = 0,
    ET_Y = 2,
    ET_X = 1,
    ET_XY = 3,
    ET_Z = 8,
    ET_ZY = 10,
    ET_ZX = 9,
    ET_ZXY = 11,
    ET_XY_XY = 7,
    ET_XY_ZXY = 15,
    ET_YZ_ZY = 42,
    ET_YZ_ZXY = 43,
    ET_YZ_XY_ZXY = 47,
    ET_ZX_ZX = 25,
    ET_ZX_ZXY = 27,
    ET_ZX_XY_ZXY = 31,
    ET_ZX_YZ_ZXY = 59,
    ET_ZX_YZ_XY_ZXY = 63
  };

  elasStiffness(daType da);
  /**
   * 	@brief		The elemental matrix-vector multiplication routine that is used
   *				by matrix-free methods. 
   * 	@param		_in	PETSc Vec which is the input vector with whom the 
   * 				product is to be calculated.
   * 	@param		_out PETSc Vec, the output of M*_in
   * 	@return		bool true if successful, false otherwise.
   *  @todo		Might have to change _in and _out to std. C arrays for speed.
   * 
   **/ 
  inline bool ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale);
  inline bool ElementalMatVec(unsigned int idx, PetscScalar *in, PetscScalar *out, double scale);

  inline bool GetElementalMatrix(int i, int j, int k, PetscScalar *mat);
  inline bool GetElementalMatrix(unsigned int idx, std::vector<ot::MatRecord>& records);

  inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale);
  inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale);

  inline bool initStencils();

  bool preMatVec();
  bool postMatVec();

  void setMu(Vec mv) {
    muVec = mv;
  }
  void setLambda(Vec lam) {
    lambdaVec = lam;
  }

  void setLame(Vec lam, Vec mu) {
    muVec = mu; lambdaVec = lam;
  }

protected:
  inline bool getElemInfo(ot::DA * da, exhaustiveElemType & sType, ot::DA::index* indices); 

private:
  void*         m_mu; 
  void*         m_lambda; 
  Vec                 muVec;     
  Vec                 lambdaVec;     

  double    m_dHx;

  double xFac, yFac, zFac;
  unsigned int maxD;

};

elasStiffness::elasStiffness(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA    = NULL;
  m_octDA   = NULL;
  m_stencil = NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();
}

bool elasStiffness::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;

  if ( m_daType == PETSC) {
    int Bijk[2][24][24] = {
      { // K1 - multiplies lambda
        {8,6,6,-8,6,6,4,-6,3,-4,-6,3,4,3,-6,-4,3,-6,2,-3,-3,-2,-3,-3},
        {6,8,6,-6,4,3,6,-8,6,-6,-4,3,3,4,-6,-3,2,-3,3,-4,-6,-3,-2,-3},
        {6,6,8,-6,3,4,3,-6,4,-3,-3,2,6,6,-8,-6,3,-4,3,-6,-4,-3,-3,-2},
        {-8,-6,-6,8,-6,-6,-4,6,-3,4,6,-3,-4,-3,6,4,-3,6,-2,3,3,2,3,3},
        {6,4,3,-6,8,6,6,-4,3,-6,-8,6,3,2,-3,-3,4,-6,3,-2,-3,-3,-4,-6},
        {6,3,4,-6,6,8,3,-3,2,-3,-6,4,6,3,-4,-6,6,-8,3,-3,-2,-3,-6,-4},
        {4,6,3,-4,6,3,8,-6,6,-8,-6,6,2,3,-3,-2,3,-3,4,-3,-6,-4,-3,-6},
        {-6,-8,-6,6,-4,-3,-6,8,-6,6,4,-3,-3,-4,6,3,-2,3,-3,4,6,3,2,3},
        {3,6,4,-3,3,2,6,-6,8,-6,-3,4,3,6,-4,-3,3,-2,6,-6,-8,-6,-3,-4},
        {-4,-6,-3,4,-6,-3,-8,6,-6,8,6,-6,-2,-3,3,2,-3,3,-4,3,6,4,3,6},
        {-6,-4,-3,6,-8,-6,-6,4,-3,6,8,-6,-3,-2,3,3,-4,6,-3,2,3,3,4,6},
        {3,3,2,-3,6,4,6,-3,4,-6,-6,8,3,3,-2,-3,6,-4,6,-3,-4,-6,-6,-8},
        {4,3,6,-4,3,6,2,-3,3,-2,-3,3,8,6,-6,-8,6,-6,4,-6,-3,-4,-6,-3},
        {3,4,6,-3,2,3,3,-4,6,-3,-2,3,6,8,-6,-6,4,-3,6,-8,-6,-6,-4,-3},
        {-6,-6,-8,6,-3,-4,-3,6,-4,3,3,-2,-6,-6,8,6,-3,4,-3,6,4,3,3,2},
        {-4,-3,-6,4,-3,-6,-2,3,-3,2,3,-3,-8,-6,6,8,-6,6,-4,6,3,4,6,3},
        {3,2,3,-3,4,6,3,-2,3,-3,-4,6,6,4,-3,-6,8,-6,6,-4,-3,-6,-8,-6},
        {-6,-3,-4,6,-6,-8,-3,3,-2,3,6,-4,-6,-3,4,6,-6,8,-3,3,2,3,6,4},
        {2,3,3,-2,3,3,4,-3,6,-4,-3,6,4,6,-3,-4,6,-3,8,-6,-6,-8,-6,-6},
        {-3,-4,-6,3,-2,-3,-3,4,-6,3,2,-3,-6,-8,6,6,-4,3,-6,8,6,6,4,3},
        {-3,-6,-4,3,-3,-2,-6,6,-8,6,3,-4,-3,-6,4,3,-3,2,-6,6,8,6,3,4},
        {-2,-3,-3,2,-3,-3,-4,3,-6,4,3,-6,-4,-6,3,4,-6,3,-8,6,6,8,6,6},
        {-3,-2,-3,3,-4,-6,-3,2,-3,3,4,-6,-6,-4,3,6,-8,6,-6,4,3,6,8,6},
        {-3,-3,-2,3,-6,-4,-6,3,-4,6,6,-8,-3,-3,2,3,-6,4,-6,3,4,6,6,8}
      },
      { // K2 - multiplies mu
        {32,6,6,-8,-6,-6,4,6,3,-10,-6,-3,4,3,6,-10,-3,-6,-4,3,3,-8,-3,-3},
        {6,32,6,6,4,3,-6,-8,-6,-6,-10,-3,3,4,6,3,-4,3,-3,-10,-6,-3,-8,-3},
        {6,6,32,6,3,4,3,6,4,3,3,-4,-6,-6,-8,-6,-3,-10,-3,-6,-10,-3,-3,-8},
        {-8,6,6,32,-6,-6,-10,6,3,4,-6,-3,-10,3,6,4,-3,-6,-8,3,3,-4,-3,-3},
        {-6,4,3,-6,32,6,6,-10,-3,6,-8,-6,-3,-4,3,-3,4,6,3,-8,-3,3,-10,-6},
        {-6,3,4,-6,6,32,-3,3,-4,-3,6,4,6,-3,-10,6,-6,-8,3,-3,-8,3,-6,-10},
        {4,-6,3,-10,6,-3,32,-6,6,-8,6,-6,-4,-3,3,-8,3,-3,4,-3,6,-10,3,-6},
        {6,-8,6,6,-10,3,-6,32,-6,-6,4,-3,3,-10,6,3,-8,3,-3,4,-6,-3,-4,-3},
        {3,-6,4,3,-3,-4,6,-6,32,6,-3,4,-3,6,-10,-3,3,-8,-6,6,-8,-6,3,-10},
        {-10,-6,3,4,6,-3,-8,-6,6,32,6,-6,-8,-3,3,-4,3,-3,-10,-3,6,4,3,-6},
        {-6,-10,3,-6,-8,6,6,4,-3,6,32,-6,-3,-8,3,-3,-10,6,3,-4,-3,3,4,-6},
        {-3,-3,-4,-3,-6,4,-6,-3,4,-6,-6,32,3,3,-8,3,6,-10,6,3,-10,6,6,-8},
        {4,3,-6,-10,-3,6,-4,3,-3,-8,-3,3,32,6,-6,-8,-6,6,4,6,-3,-10,-6,3},
        {3,4,-6,3,-4,-3,-3,-10,6,-3,-8,3,6,32,-6,6,4,-3,-6,-8,6,-6,-10,3},
        {6,6,-8,6,3,-10,3,6,-10,3,3,-8,-6,-6,32,-6,-3,4,-3,-6,4,-3,-3,-4},
        {-10,3,-6,4,-3,6,-8,3,-3,-4,-3,3,-8,6,-6,32,-6,6,-10,6,-3,4,-6,3},
        {-3,-4,-3,-3,4,-6,3,-8,3,3,-10,6,-6,4,-3,-6,32,-6,6,-10,3,6,-8,6},
        {-6,3,-10,-6,6,-8,-3,3,-8,-3,6,-10,6,-3,4,6,-6,32,3,-3,-4,3,-6,4},
        {-4,-3,-3,-8,3,3,4,-3,-6,-10,3,6,4,-6,-3,-10,6,3,32,-6,-6,-8,6,6},
        {3,-10,-6,3,-8,-3,-3,4,6,-3,-4,3,6,-8,-6,6,-10,-3,-6,32,6,-6,4,3},
        {3,-6,-10,3,-3,-8,6,-6,-8,6,-3,-10,-3,6,4,-3,3,-4,-6,6,32,-6,3,4},
        {-8,-3,-3,-4,3,3,-10,-3,-6,4,3,6,-10,-6,-3,4,6,3,-8,-6,-6,32,6,6},
        {-3,-8,-3,-3,-10,-6,3,-4,3,3,4,6,-6,-10,-3,-6,-8,-6,6,4,3,6,32,6},
        {-3,-3,-8,-3,-6,-10,-6,-3,-10,-6,-6,-8,3,3,-4,3,6,4,6,3,4,6,6,32}
      },
    };

    int*** Aijk = new int2Ptr[2];
    Aijk[0] = new int1Ptr[24];
    Aijk[1] = new int1Ptr[24];
    for (int j=0;j<24;j++) {
      Aijk[0][j] = new int[24];
      Aijk[1][j] = new int[24];
      for (int k=0;k<24;k++) {
        Aijk[0][j][k] = Bijk[0][j][k];
        Aijk[1][j][k] = Bijk[1][j][k];
      }//end k
    }//end j
    m_stencil = Aijk;

  } else {
    m_stencil = NULL;
    // allocate memory for the stencils 
    double *K = new double[2*8*18*24*24];
    // read in the stencils from the file ...
    std::ifstream in("K.inp");
    in.read((char *)K, 2*8*18*24*24*sizeof(double));
    in.close();
    m_stencil = K;
  }
  return true;
}

bool elasStiffness::preMatVec() {
  // nuVec should be set directly into matrix outside the loop ...
  if (m_daType == PETSC) {
    PetscScalar *mu; // 
    PetscScalar *lambda; // 

    CHKERRQ ( VecGetArray(muVec, &mu) );
    CHKERRQ ( VecGetArray(lambdaVec, &lambda) );
    m_mu = mu;
    m_lambda = lambda;
    // compute Hx
    PetscInt mx,my,mz;
    CHKERRQ ( DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ) ; 

    m_dHx = m_dLx/(mx -1);

    m_dHx = m_dHx; 
    m_dHx /= 72.0;
  } else {
    PetscScalar *mu; 
    PetscScalar *lambda; 
    // Get nuarray
    m_octDA->vecGetBuffer(muVec, mu, false, false, true,m_uiDof);
    m_octDA->vecGetBuffer(lambdaVec, lambda, false, false, true,m_uiDof);
    m_mu = mu;
    m_lambda = lambda;

    // compute Hx
    // For octree Hx values will change per element, so has to be 
    // computed inside the loop.
    unsigned int maxD = m_octDA->getMaxDepth();

    double xFac,yFac,zFac;

    // Get the  x,y,z factors 
    xFac = 1.0/((double)(1<<(maxD-1)));
    if (m_octDA->getDimension() > 1) {
      yFac = 1.0/((double)(1<<(maxD-1)));
      if (m_octDA->getDimension() > 2) {
        zFac = 1.0/((double)(1<<(maxD-1)));
      }
    }
  }
	return true;
}

bool elasStiffness::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  //double hy = yFac*(1<<(maxD - lev));
  //double hz = zFac*(1<<(maxD - lev));

  // @check
  double fac = hx;

  stdElemType elemType;
  ot::DA::index idx[8];

  double *K = (double *)m_stencil;

  // need child number, elemType and indices.
  alignElementAndVertices(m_octDA, elemType, idx);       

  unsigned int chNum = m_octDA->getChildNumber();
  // get A and B which are the correct 24x24 matrices.
  double *A = K + (chNum*18 + elemType)*24*24;
  double *B = K + ((8+chNum)*18 + elemType)*24*24;

  PetscScalar mu      = ((PetscScalar *) m_mu)[i];
  PetscScalar lam  = ((PetscScalar *) m_lambda)[i];

  // std::cout << "lampmu is " << lampmu <<std::endl;

  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]]   +=  fac* ( ( lam*A[24*m_uiDof*k + m_uiDof*j]   + mu*B[24*m_uiDof*k + m_uiDof*j]   )*in[m_uiDof*idx[j]] +
                                       ( lam*A[24*m_uiDof*k + m_uiDof*j+1] + mu*B[24*m_uiDof*k + m_uiDof*j+1] )*in[m_uiDof*idx[j]+1] +
                                       ( lam*A[24*m_uiDof*k + m_uiDof*j+2] + mu*B[24*m_uiDof*k + m_uiDof*j+2] )*in[m_uiDof*idx[j]+2] );

      out[m_uiDof*idx[k]+1] +=  fac* ( ( lam*A[24*(m_uiDof*k+1)+ m_uiDof*j]  + mu*B[24*(m_uiDof*k+1)+m_uiDof*j]  )*in[m_uiDof*idx[j]] +
                                       ( lam*A[24*(m_uiDof*k+1)+m_uiDof*j+1] + mu*B[24*(m_uiDof*k+1)+m_uiDof*j+1] )*in[m_uiDof*idx[j]+1] +
                                       ( lam*A[24*(m_uiDof*k+1)+m_uiDof*j+2] + mu*B[24*(m_uiDof*k+1)+m_uiDof*j+2] )*in[m_uiDof*idx[j]+2] );

      out[m_uiDof*idx[k]+2] +=  fac* ( ( lam*A[24*(m_uiDof*k+2)+m_uiDof*j]   + mu*B[24*(m_uiDof*k+2)+m_uiDof*j]   )*in[m_uiDof*idx[j]] +
                                       ( lam*A[24*(m_uiDof*k+2)+m_uiDof*j+1] + mu*B[24*(m_uiDof*k+2)+m_uiDof*j+1] )*in[m_uiDof*idx[j]+1] +
                                       ( lam*A[24*(m_uiDof*k+2)+m_uiDof*j+2] + mu*B[24*(m_uiDof*k+2)+m_uiDof*j+2] )*in[m_uiDof*idx[j]+2] );
    }//end for j
  }//end for k
  return true;
}

bool elasStiffness::GetElementalMatrix(unsigned int i, std::vector<ot::MatRecord>& records) {
	unsigned int lev = m_octDA->getLevel(i);
	double hx = xFac*(1<<(maxD - lev));
	// double hy = yFac*(1<<(maxD - lev));
	// double hz = zFac*(1<<(maxD - lev));

	PetscScalar mu      = ((PetscScalar *) m_mu)[i];
	PetscScalar lam  = ((PetscScalar *) m_lambda)[i];

	double fac = hx;

	// stdElemType elemType;
	ot::DA::index idx[8];

	unsigned int chNum = m_octDA->getChildNumber();
	m_octDA->getNodeIndices(idx); 
	unsigned char hangingMask = m_octDA->getHangingNodeIndex( m_octDA->curr() );    
	unsigned char eType = getEtype(hangingMask, chNum);

	// get A and B which are the correct 24x24 matrices.
	double *K = (double *)m_stencil;
	double *A = K + (chNum*18 + eType)*24*24;
	double *B = K + ((8+chNum)*18 + eType)*24*24;

	ot::MatRecord currRec;
	for (int k = 0;k < 8;k++) {
		for (int j=0;j<8;j++) {
			currRec.rowIdx = idx[k];
			currRec.colIdx = idx[j];
			currRec.rowDim = 0;
			currRec.colDim = 0;
			currRec.val = fac* ( ( lam*A[24*m_uiDof*k + m_uiDof*j]   +  mu*B[24*m_uiDof*k + m_uiDof*j] )) ;
			records.push_back(currRec);
			currRec.rowDim = 1;
			currRec.colDim = 1;
			records.push_back(currRec);
			currRec.rowDim = 2;
			currRec.colDim = 2;
			records.push_back(currRec);
		}//end for j
	}//end for k

	return true;
}

bool elasStiffness::GetElementalMatrix(int i, int j, int k, PetscScalar *mat){
#ifdef __DEBUG__
  std::cout << CYN"\tEntering " << __func__ << NRM << i << ", " << j << ", " << k << std::endl;
#endif

  int dof=m_uiDof;
  double fac = -m_dHx;
  int ***K = (int ***)m_stencil;
  int** A = K[0];
  int** B = K[1];

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  PetscScalar mu     = ((PetscScalar *) m_mu)[((k-z)*n + j - y)*m + i-x];
  PetscScalar lam    = ((PetscScalar *) m_lambda)[((k-z)*n + j - y)*m + i-x];

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      for (int i=0; i<dof; i++) {
        for (int j=0; j<dof; j++) {
          mat[24*(3*q+i) + 3*r+j] = fac*(lam*A[dof*q+i][dof*r+j] + mu*B[dof*q+i][dof*r+j]);
        }
      }
    }
  }
  return true;
}

bool elasStiffness::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) {
#ifdef __DEBUG__
  std::cout << CYN"\tEntering " << __func__ << NRM << i << ", " << j << ", " << k << std::endl;
#endif

  int dof=3;
  int idx[8][3]={
    {k, j, dof*i},
    {k,j,dof*(i+1)},
    {k,j+1,dof*i},
    {k,j+1,dof*(i+1)},
    {k+1, j, dof*i},
    {k+1,j,dof*(i+1)},
    {k+1,j+1,dof*i},
    {k+1,j+1,dof*(i+1)}               
  };             

  double fac = -m_dHx*scale;
  int ***K = (int ***)m_stencil;
  int** A = K[0];
  int** B = K[1];

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  PetscScalar mu     = ((PetscScalar *) m_mu)[((k-z)*n + j - y)*m + i-x];
  PetscScalar lam    = ((PetscScalar *) m_lambda)[((k-z)*n + j - y)*m + i-x];

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      out[idx[q][0]][idx[q][1]][idx[q][2]]   +=  fac* ( ( lam*A[dof*q][dof*r] + mu*B[dof*q][dof*r] )*in[idx[r][0]][idx[r][1]][idx[r][2]] +
                                                        ( lam*A[dof*q][dof*r+1] + mu*B[dof*q][dof*r+1] )*in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                        ( lam*A[dof*q][dof*r+2] + mu*B[dof*q][dof*r+2] )*in[idx[r][0]][idx[r][1]][idx[r][2]+2] );

      out[idx[q][0]][idx[q][1]][idx[q][2]+1] +=  fac* ( ( lam*A[dof*q+1][dof*r]   + mu*B[dof*q+1][dof*r]    ) * in[idx[r][0]][idx[r][1]][idx[r][2]] +
                                                        ( lam*A[dof*q+1][dof*r+1] + mu*B[dof*q+1][dof*r+1]  ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                        ( lam*A[dof*q+1][dof*r+2] + mu*B[dof*q+1][dof*r+2]  ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );

      out[idx[q][0]][idx[q][1]][idx[q][2]+2] +=   fac* ( ( lam*A[dof*q+2][dof*r]   + mu*B[dof*q+2][dof*r]   ) * in[idx[r][0]][idx[r][1]][idx[r][2]] +
                                                         ( lam*A[dof*q+2][dof*r+1] + mu*B[dof*q+2][dof*r+1] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                         ( lam*A[dof*q+2][dof*r+2] + mu*B[dof*q+2][dof*r+2] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );
      /*
      out[idx[q][0]][idx[q][1]][idx[q][2]]   +=  fac* ( ( mu*A[dof*q][dof*r] + lampmu*B[dof*q][dof*r] )*in[idx[r][0]][idx[r][1]][idx[r][2]] +
                                                        ( mu*A[dof*q][dof*r+1] + lampmu*B[dof*q][dof*r+1] )*in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                        ( mu*A[dof*q][dof*r+2] + lampmu*B[dof*q][dof*r+2] )*in[idx[r][0]][idx[r][1]][idx[r][2]+2] );

      out[idx[q][0]][idx[q][1]][idx[q][2]+1] +=  fac* ( ( mu*A[dof*q+1][dof*r] + lampmu*B[dof*q+1][dof*r] )*(in[idx[r][0]][idx[r][1]][idx[r][2]]) +
                                                        ( mu*A[dof*q+1][dof*r+1] + lampmu*B[dof*q+1][dof*r+1] )*in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                        ( mu*A[dof*q+1][dof*r+2] + lampmu*B[dof*q+1][dof*r+2] )*in[idx[r][0]][idx[r][1]][idx[r][2]+2] );

      out[idx[q][0]][idx[q][1]][idx[q][2]+2]   +=  fac* ( ( mu*A[dof*q+2][dof*r] + lampmu*B[dof*q+2][dof*r] )*in[idx[r][0]][idx[r][1]][idx[r][2]] +
                                                          ( mu*A[dof*q+2][dof*r+1] + lampmu*B[dof*q+2][dof*r+1] )*in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                          ( mu*A[dof*q+2][dof*r+2] + lampmu*B[dof*q+2][dof*r+2] )*in[idx[r][0]][idx[r][1]][idx[r][2]+2] );
      */
    }
  }
}

bool elasStiffness::postMatVec() {
  if ( m_daType == PETSC) {
    PetscScalar *mu = (PetscScalar *)m_mu;
    CHKERRQ( VecRestoreArray(muVec, &mu) );
    PetscScalar *lam = (PetscScalar *)m_lambda;
    CHKERRQ( VecRestoreArray(lambdaVec, &lam) );
  } else {
    PetscScalar *mu = (PetscScalar *)m_mu;
    m_octDA->vecRestoreBuffer(muVec, mu, false, false, true, m_uiDof);
    PetscScalar *lam = (PetscScalar *)m_lambda;
    m_octDA->vecRestoreBuffer(lambdaVec, lam, false, false, true, m_uiDof);
  }
  return true;
}

bool elasStiffness::getElemInfo(ot::DA* da, exhaustiveElemType& sType, ot::DA::index *indices) {
  sType = ET_N;
  da->getNodeIndices(indices); 

  if (da->isHanging(da->curr())) {

    int childNum = da->getChildNumber();
    Point pt = da->getCurrentOffset();   

    unsigned char hangingMask = da->getHangingNodeIndex(da->curr());    

    //Change HangingMask and indices based on childNum
    mapVtxAndFlagsToOrientation(childNum, indices, hangingMask);    

    unsigned char eType = ((126 & hangingMask)>>1);

    reOrderIndices(eType, indices);
  }//end if hangingElem.
  return true;
}

bool elasStiffness::ElementalMatGetDiagonal(unsigned int i, PetscScalar *diag, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  // double hy = yFac*(1<<(maxD - lev));
  // double hz = zFac*(1<<(maxD - lev));

  double fac = hx;

  stdElemType elemType;
  ot::DA::index idx[8];

  float *K = (float *)m_stencil;

  // need child number, elemType and indices.
  alignElementAndVertices(m_octDA, elemType, idx);       

  unsigned int chNum = m_octDA->getChildNumber();
  // get A and B which are the correct 24x24 matrices.
  float *A = K + (chNum*18 + elemType)*24*24;
  float *B = K + ((8+chNum)*18 + elemType)*24*24;

  PetscScalar mu      = ((PetscScalar *) m_mu)[i];
  PetscScalar lam  = ((PetscScalar *) m_lambda)[i];

  // std::cout << "lampmu is " << lampmu <<std::endl;

  for (int k = 0;k < 8;k++) {
    diag[m_uiDof*idx[k]]   +=  fac* ( lam*A[24*m_uiDof*k + m_uiDof*k]   + mu*B[24*m_uiDof*k + m_uiDof*k]   ) ;

    diag[m_uiDof*idx[k]+1] +=  fac* ( lam*A[24*(m_uiDof*k+1)+m_uiDof*k+1] + mu*B[24*(m_uiDof*k+1)+m_uiDof*k+1] ) ;

    diag[m_uiDof*idx[k]+2] +=  fac* ( lam*A[24*(m_uiDof*k+2)+m_uiDof*k+2] + mu*B[24*(m_uiDof*k+2)+m_uiDof*k+2] ) ;
  }//end for k
  return true;
}

bool elasStiffness::ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) {
#ifdef __DEBUG__
  std::cout << CYN"\tEntering " << __func__ << NRM << i << ", " << j << ", " << k << std::endl;
#endif

  int dof=3;
  int idx[8][3]={
    {k, j, dof*i},
    {k,j,dof*(i+1)},
    {k,j+1,dof*i},
    {k,j+1,dof*(i+1)},
    {k+1, j, dof*i},
    {k+1,j,dof*(i+1)},
    {k+1,j+1,dof*i},
    {k+1,j+1,dof*(i+1)}               
  };             

  double fac = -m_dHx*scale;
  int ***K = (int ***)m_stencil;
  int** A = K[0];
  int** B = K[1];

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );
  // std::cout << "fac is " << fac << std::endl;

  PetscScalar mu     = ((PetscScalar *) m_mu)[((k-z)*n + j - y)*m + i-x];
  PetscScalar lam    = ((PetscScalar *) m_lambda)[((k-z)*n + j - y)*m + i-x];
  // std::cout << "Mu and Lambda are " << mu << ", " << lam << std::endl; 

  for (int q = 0; q < 8; q++) {

    diag[idx[q][0]][idx[q][1]][idx[q][2]]   +=  fac* ( lam*A[dof*q][dof*q] + mu*B[dof*q][dof*q] );
    diag[idx[q][0]][idx[q][1]][idx[q][2]+1] +=  fac* ( lam*A[dof*q+1][dof*q+1] + mu*B[dof*q+1][dof*q+1]  ) ;
    diag[idx[q][0]][idx[q][1]][idx[q][2]+2] +=  fac* ( lam*A[dof*q+2][dof*q+2] + mu*B[dof*q+2][dof*q+2]  );
    /*
    diag[idx[q][0]][idx[q][1]][idx[q][2]]   += fac* ( mu*A[dof*q][dof*q] + lampmu*B[dof*q][dof*q] ) ;
    diag[idx[q][0]][idx[q][1]][idx[q][2]+1] +=  fac* ( mu*A[dof*q+1][dof*q+1] + lampmu*B[dof*q+1][dof*q+1] ) ;
    diag[idx[q][0]][idx[q][1]][idx[q][2]+2] +=  fac* ( mu*A[dof*q+2][dof*q+2] + lampmu*B[dof*q+2][dof*q+2] ) ;
    */

  }
}

#endif /*_STIFFNESSMATRIX_H_*/

