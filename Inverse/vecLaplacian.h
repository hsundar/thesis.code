/**
 *  @file	vecLaplacian.h
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
#ifndef _VEC_LAPLACIAN_H_
#define _VEC_LAPLACIAN_H_

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

class vecLaplacian : public feMatrix<vecLaplacian> {
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

  vecLaplacian(daType da);
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

  inline bool initStencils();

  bool preMatVec();
  bool postMatVec();

protected:
  inline bool getElemInfo(ot::DA * da, exhaustiveElemType & sType, ot::DA::index* indices); 

private:

  double    m_dHx;

  double xFac, yFac, zFac;
  unsigned int maxD;
};

vecLaplacian::vecLaplacian(daType da) {
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

bool vecLaplacian::initStencils() {
  typedef int* int1Ptr;
  typedef int** int2Ptr;

  if ( m_daType == PETSC) {

    int Bjk[8][8] = {
      { 64,   0,   0, -16,   0, -16, -16, -16},
      {0,  64, -16,   0, -16,   0, -16, -16},
      {0, -16,  64,   0, -16, -16,   0, -16},
      { -16,   0,   0,  64, -16, -16, -16,   0},
      { 0, -16, -16, -16,  64,   0,   0, -16},
      { -16,   0, -16, -16,   0,  64, -16,   0},
      { -16, -16,   0, -16,   0, -16,  64,   0},
      { -16, -16, -16,   0, -16,   0,   0,  64} 
      /*
      {24,0,0,0,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6,0,0,-6,0,0},
      {0,24,0,0,0,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6,0,0,-6,0},
      {0,0,24,0,0,0,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6,0,0,-6},
      {0,0,0,24,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6,0,0},
      {0,0,0,0,24,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6,0},
      {0,0,0,0,0,24,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,-6},
      {0,0,0,-6,0,0,24,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,-6,0,0},
      {0,0,0,0,-6,0,0,24,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,-6,0},
      {0,0,0,0,0,-6,0,0,24,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,-6},
      {-6,0,0,0,0,0,0,0,0,24,0,0,-6,0,0,-6,0,0,-6,0,0,0,0,0},
      {0,-6,0,0,0,0,0,0,0,0,24,0,0,-6,0,0,-6,0,0,-6,0,0,0,0},
      {0,0,-6,0,0,0,0,0,0,0,0,24,0,0,-6,0,0,-6,0,0,-6,0,0,0},
      {0,0,0,-6,0,0,-6,0,0,-6,0,0,24,0,0,0,0,0,0,0,0,-6,0,0},
      {0,0,0,0,-6,0,0,-6,0,0,-6,0,0,24,0,0,0,0,0,0,0,0,-6,0},
      {0,0,0,0,0,-6,0,0,-6,0,0,-6,0,0,24,0,0,0,0,0,0,0,0,-6},
      {-6,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,24,0,0,-6,0,0,0,0,0},
      {0,-6,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,24,0,0,-6,0,0,0,0},
      {0,0,-6,0,0,0,0,0,-6,0,0,-6,0,0,0,0,0,24,0,0,-6,0,0,0},
      {-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,24,0,0,0,0,0},
      {0,-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,24,0,0,0,0},
      {0,0,-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,-6,0,0,24,0,0,0},
      {-6,0,0,-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,0,0,0,24,0,0},
      {0,-6,0,0,-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,0,0,0,24,0},
      {0,0,-6,0,0,-6,0,0,-6,0,0,0,0,0,-6,0,0,0,0,0,0,0,0,24}
      */
    };

    int** Ajk = new int1Ptr[8];
    for (int j=0;j<8;j++) {
      Ajk[j] = new int[8];
      for (int k=0;k<8;k++) {
        Ajk[j][k] = Bjk[j][k];
      }//end k
    }//end j
    m_stencil = Ajk;

    /*
    int Bijk[24][24] = {
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
    };
   

    int** Aijk = new int1Ptr[24];
    for (int j=0;j<24;j++) {
      Aijk[j] = new int[24];
      for (int k=0;k<24;k++) {
        Aijk[j][k] = Bijk[j][k];
      }//end k
    }//end j
    m_stencil = Aijk;
    */

  } else {
    m_stencil = NULL;
    // allocate memory for the stencils 
    float *K = new float[8*18*24*24];
    // read in the stencils from the file ...
    std::ifstream in("vLap.dat");
    in.read((char *)K, 8*18*24*24);
    in.close();
    m_stencil = K;
  }
  return true;
}

bool vecLaplacian::preMatVec() {
  // nuVec should be set directly into matrix outside the loop ...
  if (m_daType == PETSC) {
    // compute Hx
    PetscInt mx,my,mz;
    CHKERRQ ( DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0) ) ; 

    m_dHx = m_dLx/(mx -1);

    m_dHx /= 72.0;
  } else {
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
}

bool vecLaplacian::ElementalMatVec(unsigned int i, PetscScalar *in, PetscScalar *out, double scale) {
  unsigned int lev = m_octDA->getLevel(i);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  // @check
  double fac = -hx*scale;

  stdElemType elemType;
  ot::DA::index idx[8];

  float *K = (float *)m_stencil;

  // need child number, elemType and indices.
  alignElementAndVertices(m_octDA, elemType, idx);       

  unsigned int chNum = m_octDA->getChildNumber();
  // get A and B which are the correct 24x24 matrices.
  float *A = K + (chNum*18 + elemType)*24*24;

  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      out[m_uiDof*idx[k]]   +=  fac* ( ( A[24*m_uiDof*k + m_uiDof*j]   )*in[m_uiDof*idx[j]] +
                                       ( A[24*m_uiDof*k + m_uiDof*j+1] )*in[m_uiDof*idx[j]+1] +
                                       ( A[24*m_uiDof*k + m_uiDof*j+2] )*in[m_uiDof*idx[j]+2] );

      out[m_uiDof*idx[k]+1] +=  fac* ( ( A[24*(m_uiDof*k+1)+ m_uiDof*j]   )*in[m_uiDof*idx[j]] +
                                       ( A[24*(m_uiDof*k+1)+m_uiDof*j+1]  )*in[m_uiDof*idx[j]+1] +
                                       ( A[24*(m_uiDof*k+1)+m_uiDof*j+2]  )*in[m_uiDof*idx[j]+2] );

      out[m_uiDof*idx[k]+2] +=  fac* ( ( A[24*(m_uiDof*k+2)+m_uiDof*j]  )*in[m_uiDof*idx[j]] +
                                       ( A[24*(m_uiDof*k+2)+m_uiDof*j+1])*in[m_uiDof*idx[j]+1] +
                                       ( A[24*(m_uiDof*k+2)+m_uiDof*j+2])*in[m_uiDof*idx[j]+2] );
    }//end for j
  }//end for k
  return true;
}

bool vecLaplacian::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale){
#ifdef __DEBUG__  
  std::cout << "Entering " << __func__ <<  " " << i << ", " << j << ", " << k << std::endl;
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

  int** A = (int **)m_stencil;

  int x,y,z,m,n,p;
  CHKERRQ( DAGetCorners(m_DA, &x, &y, &z, &m, &n, &p) );

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      out[idx[q][0]][idx[q][1]][idx[q][2]]     +=  fac* ( ( A[q][r]     ) * in[idx[r][0]][idx[r][1]][idx[r][2]] ); // +
                                                       //  ( A[dof*q][dof*r+1]   ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                       //   ( A[dof*q][dof*r+2]   ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );

      /*
      out[idx[q][0]][idx[q][1]][idx[q][2]+1]   +=  fac* ( ( A[dof*q+1][dof*r]   ) * in[idx[r][0]][idx[r][1]][idx[r][2]]   +
                                                          ( A[dof*q+1][dof*r+1] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                          ( A[dof*q+1][dof*r+2] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );
      out[idx[q][0]][idx[q][1]][idx[q][2]+2]   +=  fac* ( ( A[dof*q+2][dof*r]   ) * in[idx[r][0]][idx[r][1]][idx[r][2]]   +
                                                          ( A[dof*q+2][dof*r+1] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] +
                                                          ( A[dof*q+2][dof*r+2] ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );
      *
      out[idx[q][0]][idx[q][1]][idx[q][2]]     +=  fac* ( ( A[q][r]  ) * in[idx[r][0]][idx[r][1]][idx[r][2]] );
      out[idx[q][0]][idx[q][1]][idx[q][2]+1]     +=  fac* ( ( A[q][r]  ) * in[idx[r][0]][idx[r][1]][idx[r][2]+1] );
      out[idx[q][0]][idx[q][1]][idx[q][2]+2]     +=  fac* ( ( A[q][r]  ) * in[idx[r][0]][idx[r][1]][idx[r][2]+2] );
      */
    }
  }
}

bool vecLaplacian::postMatVec() {
  if ( m_daType == PETSC) {

  } else {

  }
  return true;
}

bool vecLaplacian::getElemInfo(ot::DA* da, exhaustiveElemType& sType, ot::DA::index *indices) {
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

#endif /*_STIFFNESSMATRIX_H_*/

