
/**
*  @file	elasMass.h
*  @brief	Main class for finite element assembly of a generic
*           mass matrix.
*  @author	Hari Sundar
*  @date	5/8/07
* 
*  Main class for finite element assembly of a generic mass matrix. The density 
*  image needs to be set and is used for the computation of the matrix. MatVec 
*  functions are virtual and in most cases need not be specified as they are 
*  obtained from the parent feMatrix class. In most cases only the assembly of 
*  the element matrix needs to be done within this or a derived class.
*/

#ifndef _MASSMATRIX_H_
#define _MASSMATRIX_H_

/**
*  @brief	Main class for finite element assembly of a generic mass matrix.
*  @author	Hari Sundar
*  @date	5/8/07
*
*  Main class for finite element assembly of a generic mass matrix. The density 
*  image needs to be set and is used for the computation of the matrix. MatVec 
*  functions are virtual and in most cases need not be specified as they are 
*  obtained from the parent feMatrix class. In most cases only the assembly of 
*  the element matrix needs to be done within this or a derived class.
*/

#include "feMatrix.h"

class elasMass : public feMatrix<elasMass> {
public:
  elasMass(daType da);
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
  inline bool GetElementalMatrix(unsigned int idx, std::vector<ot::MatRecord> &record);

  inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale);
  inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale);

  inline bool initStencils();

  bool preMatVec();
  bool postMatVec();

  void setDensity(Vec rho) {
    rhoVec = rho;
  }

private:
  // Density vector
  void *          m_rho; 
  Vec             rhoVec;     

  double          m_dHx;

  double xFac, yFac, zFac;
  unsigned int maxD;
};


elasMass::elasMass(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA    = NULL;
  m_octDA   = NULL;
  m_stencil = NULL;
  m_rho = NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();
}

bool elasMass::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;

  if (m_daType == PETSC) {
    int Bjk[8][8] =    {
      { 64, 32, 32, 16, 32, 16, 16,  8},
      { 32, 64, 16, 32, 16, 32,  8, 16},
      { 32, 16, 64, 32, 16,  8, 32, 16},
      { 16, 32, 32, 64,  8, 16, 16, 32},
      { 32, 16, 16,  8, 64, 32, 32, 16},
      { 16, 32,  8, 16, 32, 64, 16, 32},
      { 16,  8, 32, 16, 32, 16, 64, 32},
      {  8, 16, 16, 32, 16, 32, 32, 64}
    };

    int** Ajk = new int1Ptr[8];
    for (int j=0;j<8;j++) {
      Ajk[j] = new int[8];
      for (int k=0;k<8;k++) {
        Ajk[j][k] = Bjk[j][k];
      }//end k
    }//end j
    m_stencil = Ajk;
  } else {
    int Bijk[8][8][8] = {
      //Type-0:No Hanging
      {{ 64, 32, 32, 16, 32, 16, 16,  8},
        { 32, 64, 16, 32, 16, 32,  8, 16},
        { 32, 16, 64, 32, 16,  8, 32, 16},
        { 16, 32, 32, 64,  8, 16, 16, 32},
        { 32, 16, 16,  8, 64, 32, 32, 16},
        { 16, 32,  8, 16, 32, 64, 16, 32},
        { 16,  8, 32, 16, 32, 16, 64, 32},
        {  8, 16, 16, 32, 16, 32, 32, 64}},

      //Type-1: Y Hanging
      {{ 112,  40,  32,  32,  40,  20,  32,  16},
        {  40,  64,   8,  32,  16,  32,   8,  16},
        {  32,   8,  16,  16,   8,   4,  16,   8},
        {  32,  32,  16,  64,   8,  16,  16,  32},
        {  40,  16,   8,   8,  64,  32,  32,  16},
        {  20,  32,   4,  16,  32,  64,  16,  32},
        {  32,   8,  16,  16,  32,  16,  64,  32},
        {  16,  16,   8,  32,  16,  32,  32,  64}},

      //Type-2: X and Y Hanging
      {{ 168,  36,  36,  48,  48,  36,  36,  24},
        {  36,  16,   4,  16,   8,  16,   4,   8},
        {  36,   4,  16,  16,   8,   4,  16,   8},
        {  48,  16,  16,  64,   8,  16,  16,  32},
        {  48,   8,   8,   8,  64,  32,  32,  16},
        {  36,  16,   4,  16,  32,  64,  16,  32},
        {  36,   4,  16,  16,  32,  16,  64,  32},
        {  24,   8,   8,  32,  16,  32,  32,  64}},

      //Type-3: X and Y and Z Hanging
      {{ 232,  40,  40,  52,  40,  52,  52,  32},
        {  40,  16,   4,  16,   4,  16,   4,   8},
        {  40,   4,  16,  16,   4,   4,  16,   8},
        {  52,  16,  16,  64,   4,  16,  16,  32},
        {  40,   4,   4,   4,  16,  16,  16,   8},
        {  52,  16,   4,  16,  16,  64,  16,  32},
        {  52,   4,  16,  16,  16,  16,  64,  32},
        {  32,   8,   8,  32,   8,  32,  32,  64}},

      //Type-4:XY and X and Y Hanging
      {{ 196,  56,  56,  16,  50,  40,  40,  32},
        {  56,  28,  16,   8,  10,  20,   8,  16},
        {  56,  16,  28,   8,  10,   8,  20,  16},
        {  16,   8,   8,   4,   2,   4,   4,   8},
        {  50,  10,  10,   2,  64,  32,  32,  16},
        {  40,  20,   8,   4,  32,  64,  16,  32},
        {  40,   8,  20,   4,  32,  16,  64,  32},
        {  32,  16,  16,   8,  16,  32,  32,  64}},

      //Type-5:XY and X and Y and Z Hanging
      {{ 262,  61,  61,  17,  41,  56,  56,  40},
        {  61,  28,  16,   8,   5,  20,   8,  16},
        {  61,  16,  28,   8,   5,   8,  20,  16},
        {  17,   8,   8,   4,   1,   4,   4,   8},
        {  41,   5,   5,   1,  16,  16,  16,   8},
        {  56,  20,   8,   4,  16,  64,  16,  32},
        {  56,   8,  20,   4,  16,  16,  64,  32},
        {  40,  16,  16,   8,   8,  32,  32,  64}},

      //Type-6:XY and YZ and X and Y and Z Hanging
      {{ 294,  63,  84,  18,  63,  60,  18,  48},
        {  63,  28,  18,   8,   7,  20,   2,  16},
        {  84,  18,  42,   9,  18,  12,   9,  24},
        {  18,   8,   9,   4,   2,   4,   1,   8},
        {  63,   7,  18,   2,  28,  20,   8,  16},
        {  60,  20,  12,   4,  20,  64,   4,  32},
        {  18,   2,   9,   1,   8,   4,   4,   8},
        {  48,  16,  24,   8,  16,  32,   8,  64}},

      //Type-7: All 6 Hanging
      {{ 328,  87,  87,  19,  87,  19,  19,  56},
        {  87,  42,  21,   9,  21,   9,   3,  24},
        {  87,  21,  42,   9,  21,   3,   9,  24},
        {  19,   9,   9,   4,   3,   1,   1,   8},
        {  87,  21,  21,   3,  42,   9,   9,  24},
        {  19,   9,   3,   1,   9,   4,   1,   8},
        {  19,   3,   9,   1,   9,   1,   4,   8},
        {  56,  24,  24,   8,  24,   8,   8,  64}}
    };

    int ***Aijk = new int2Ptr[8];
    for (int i=0;i<8;i++) {
      Aijk[i] = new int1Ptr[8];
      for (int j=0;j<8;j++) {
        Aijk[i][j] = new int[8];
        for (int k=0;k<8;k++) {
          Aijk[i][j][k] = Bijk[i][j][k];
        }//end k
      }//end j
    }//end i
    m_stencil = Aijk;
  }
  return true;
}

bool elasMass::preMatVec() {
  // std::cout << "Entering " << __func__ << std::endl;
  if (m_daType == PETSC) {
    // compute Hx
    PetscInt mx,my,mz;

    PetscScalar *rho; // 
    CHKERRQ( VecGetArray(rhoVec, &rho) );
    m_rho = rho;

    CHKERRQ( DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ); 
    m_dHx = m_dLx/(mx -1);
    m_dHx = m_dHx*m_dHx*m_dHx;
    m_dHx /= 1728.0;

    // std::cout << "Hx is " << m_dHx << std::endl;
  } else {
    unsigned int maxD = m_octDA->getMaxDepth();

    double xFac,yFac,zFac;

    PetscScalar *rho; 
    // Get nuarray
    m_octDA->vecGetBuffer(rhoVec, rho, false, false, true, m_uiDof);
    m_rho = rho;

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
  // std::cout << "Leaving " << __func__ << std::endl;
}

bool elasMass::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  PetscScalar *rho  = (PetscScalar *) m_rho;

  double fac = rho[i]*scale*hx*hy*hz/1728.0;

  stdElemType elemType;
  ot::DA::index idx[8];

  int ***Aijk = (int ***)m_stencil;

  alignElementAndVertices(m_octDA, elemType, idx);       


  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]]   += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]];
      out[m_uiDof*idx[k]+1] += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]+1];
      out[m_uiDof*idx[k]+2] += fac*(Aijk[elemType][k][j])*in[m_uiDof*idx[j]+2];
    }//end for j
  }//end for k

  return true;
}

bool elasMass::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale){
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

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  PetscScalar *rho  = (PetscScalar *) m_rho;
  double stencilScale =  m_dHx*scale*rho[((k-z)*n + j - y)*m + i-x];
  int **Ajk = (int **)m_stencil;

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) { 
      out[idx[q][0]][idx[q][1]][idx[q][2]]   += stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]];
      out[idx[q][0]][idx[q][1]][idx[q][2]+1] += stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]+1];
      out[idx[q][0]][idx[q][1]][idx[q][2]+2] += stencilScale*Ajk[q][r]*in[idx[r][0]][idx[r][1]][idx[r][2]+2];
    }
  }
}

  bool elasMass::GetElementalMatrix(unsigned int i, std::vector<ot::MatRecord> &records) {
		unsigned int lev = m_octDA->getLevel(i);
	  double hx = xFac*(1<<(maxD - lev));
	  double hy = yFac*(1<<(maxD - lev));
	  double hz = zFac*(1<<(maxD - lev));

	  PetscScalar *rho  = (PetscScalar *) m_rho;

	  double fac = rho[i]*hx*hy*hz/1728.0;

	  stdElemType elemType;
	  ot::DA::index idx[8];

	  int ***Aijk = (int ***)m_stencil;

	  alignElementAndVertices(m_octDA, elemType, idx);       

		ot::MatRecord currRec;
	  for (int k = 0;k < 8;k++) {
	    for (int j=0;j<8;j++) {
				currRec.rowIdx = idx[k];
				currRec.colIdx = idx[j];
				currRec.rowDim = 0;
				currRec.colDim = 0;
				currRec.val = fac*(Aijk[elemType][k][j]);
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

bool elasMass::GetElementalMatrix(int i, int j, int k, PetscScalar *mat){
  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  PetscScalar *rho  = (PetscScalar *) m_rho;
  double stencilScale =  m_dHx*rho[((k-z)*n + j - y)*m + i-x];
  int **Ajk = (int **)m_stencil;

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      mat[24*3*q + 3*r] = stencilScale*Ajk[q][r];
      mat[24*(3*q+1) + 3*r+1] = stencilScale*Ajk[q][r];
      mat[24*(3*q+2) + 3*r+2] = stencilScale*Ajk[q][r];
    }
  }
  return true;
}

bool elasMass::postMatVec() {
  // std::cout << "Entering " << __func__ << std::endl;
  if ( m_daType == PETSC) {
    PetscScalar *rho = (PetscScalar *)m_rho;
    int ierr = VecRestoreArray(rhoVec, &rho);
    CHKERRQ(ierr);
  } else {
    PetscScalar *rho = (PetscScalar *)m_rho;
    m_octDA->vecRestoreBuffer(rhoVec, rho, false, false, true,m_uiDof);
  }
  // std::cout << "Leaving " << __func__ << std::endl;
  return true;
}

bool elasMass::ElementalMatGetDiagonal(unsigned int i, PetscScalar *diag, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  PetscScalar *rho  = (PetscScalar *) m_rho;

  double fac = rho[i]*scale*hx*hy*hz/1728.0;

  stdElemType elemType;
  ot::DA::index idx[8];

  int ***Aijk = (int ***)m_stencil;

  alignElementAndVertices(m_octDA, elemType, idx);       


  for (int k = 0;k < 8;k++) {
      diag[m_uiDof*idx[k]]   += fac*(Aijk[elemType][k][k]);
      diag[m_uiDof*idx[k]+1] += fac*(Aijk[elemType][k][k]);
      diag[m_uiDof*idx[k]+2] += fac*(Aijk[elemType][k][k]);
  }//end for k

  return true;
}

bool elasMass::ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale){
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

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  PetscScalar *rho  = (PetscScalar *) m_rho;
  double stencilScale =  m_dHx*scale*rho[((k-z)*n + j - y)*m + i-x];
  int **Ajk = (int **)m_stencil;

  for (int q = 0; q < 8; q++) {
      diag[idx[q][0]][idx[q][1]][idx[q][2]]   += stencilScale*Ajk[q][q];
      diag[idx[q][0]][idx[q][1]][idx[q][2]+1] += stencilScale*Ajk[q][q];
      diag[idx[q][0]][idx[q][1]][idx[q][2]+2] += stencilScale*Ajk[q][q];
  }
}

#endif /*_MASSMATRIX_H_*/

