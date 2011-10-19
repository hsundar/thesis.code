/**
 *  @file	cardiacDynamic.h
 *  @brief	Main class to set the activation force within the
 *         myocardium.
 *  @author Hari Sundar
 *  @date	  5/15/7
 * 
 *  Main  class to set  the activation force within the
 *  myocardium.
 */

#ifndef __CARDIAC_DYNAMIC_H_
#define __CARDIAC_DYNAMIC_H_

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
class cardiacDynamic : public feVector<cardiacDynamic> {
public: 
  cardiacDynamic(daType da);

  bool initStencils();

  inline bool ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale);
  inline bool ElementalAddVec(unsigned int index, PetscScalar *in, double scale);

  bool preAddVec();
  bool postAddVec();

  void setFDynamic(std::vector<Vec> nv) {
      m_fdynamic = nv;
    }

  private:
	 void*    m_farray;
	 std::vector<Vec>      m_fdynamic;
	 Vec       m_vecFLocal;

     double      xFac,yFac,zFac;
	 unsigned int maxD;
     double 		 m_dHx;
};

cardiacDynamic::cardiacDynamic(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA    = NULL;
  m_octDA   = NULL;
  m_stencil = NULL;

  m_farray = NULL;

  // initialize the stencils ...
  initStencils();
  if (da == OCT)
    initOctLut();

}
bool cardiacDynamic::initStencils() {
#ifdef __DEBUG__    
  std::cout << "Entering " << __func__ << std::endl;
#endif    
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
#ifdef __DEBUG__  
  std::cout << "Leaving " << __func__ << std::endl;
#endif  
  return true;
}

bool cardiacDynamic::preAddVec() {
#ifdef __DEBUG__    
  std::cout << "Entering " << __func__ << std::endl;
#endif    
  if (m_daType == PETSC){
  int ierr;
  // compute Hx
  PetscInt mx,my,mz;

  ierr = DAGetInfo(m_DA,0, &mx, &my, &mz, 0,0,0,0,0,0,0); CHKERRQ(ierr);

  m_dHx = m_dLx/(mx -1);

  m_dHx = m_dHx*m_dHx*m_dHx;

  PetscScalar ***farray;
  // double norm;
#ifdef __DEBUG__
  VecNorm(m_fdynamic[m_iCurrentDynamicIndex], NORM_INFINITY, &norm);
  std::cout << " norm of dynamic @ step " << m_time->currentstep << " = " << norm << std::endl;
#endif

  // std::cout << "Current Index is " << m_iCurrentDynamicIndex << std::endl;
  ierr = DAGetLocalVector(m_DA,&m_vecFLocal); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(m_DA,m_fdynamic[m_iCurrentDynamicIndex],INSERT_VALUES,m_vecFLocal); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(m_DA,m_fdynamic[m_iCurrentDynamicIndex],INSERT_VALUES,m_vecFLocal); CHKERRQ(ierr);

  //  ierr = DAVecGetArray(m_DA,m_fdynamic[m_time->currentstep - 1],&farray);
  ierr = DAVecGetArray(m_DA,m_vecFLocal,&farray);
  m_farray = farray;
  
  m_dHx /= 1728.0;
  // m_dHx = 1.0/1728.0;
  
  CHKERRQ(ierr);
  }else{
	 unsigned int maxD = m_octDA->getMaxDepth();

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

bool cardiacDynamic::postAddVec() {
#ifdef __DEBUG__
  std::cout << "Entering " << __func__ << std::endl;
#endif
  if ( m_daType == PETSC) {
    PetscScalar ***farray = (PetscScalar ***)m_farray;
    int ierr;
    //  ierr = DAVecRestoreArray(m_DA,m_fdynamic[m_time->currentstep - 1],&farray); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(m_DA,m_vecFLocal,&farray); CHKERRQ(ierr);
    ierr = DARestoreLocalVector(m_DA,&m_vecFLocal); CHKERRQ(ierr);
  } else {
    
  }
#ifdef __DEBUG__  
  std::cout << "Leaving " << __func__ << std::endl;
#endif  
  return true;
}

bool cardiacDynamic::ElementalAddVec(int i, int j, int k, PetscScalar ***in, double scale) {
#ifdef __DEBUG__
    std::cout << CYN"\tEntering " << __func__ << NRM << i <<", " << j << ", " << k << std::endl;
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

  PetscScalar ***fibers = (PetscScalar ***)m_farray;

  double stencilScale =  m_dHx*scale;
  int **Ajk = (int **)m_stencil;

  for (int q = 0; q < 8; q++) {
    for (int r = 0; r < 8; r++) {
      in[idx[q][0]][idx[q][1]][idx[q][2]]   += stencilScale*Ajk[q][r] * fibers[idx[r][0]][idx[r][1]][idx[r][2]];  
      in[idx[q][0]][idx[q][1]][idx[q][2]+1] += stencilScale*Ajk[q][r] * fibers[idx[r][0]][idx[r][1]][idx[r][2]+1];
      in[idx[q][0]][idx[q][1]][idx[q][2]+2] += stencilScale*Ajk[q][r] * fibers[idx[r][0]][idx[r][1]][idx[r][2]+2];

      // in[idx[q][0]][idx[q][1]][idx[q][2]]   += scale* fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]];  
      // in[idx[q][0]][idx[q][1]][idx[q][2]+1] += scale* fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+1];
      // in[idx[q][0]][idx[q][1]][idx[q][2]+2] += scale* fibers[idx[q][0]][idx[q][1]][dof*idx[q][2]+2];
    }
  }

  //std::cout << "Stencil scale is " << stencilScale << std::endl;
#ifdef __DEBUG__    
//  std::cout << "Leaving " << __func__ << std::endl;
#endif    
  return true;
}

bool cardiacDynamic::ElementalAddVec(unsigned int index, PetscScalar *in, double scale) {
#ifdef __DEBUG__  
  std::cout << "Entering " << __func__ << std::endl;
#endif  
  unsigned int lev = m_octDA->getLevel(index);
  double hx = xFac*(1<<(maxD - lev));
  double hy = yFac*(1<<(maxD - lev));
  double hz = zFac*(1<<(maxD - lev));

  double fac = scale*hx*hy*hz/1728.0;

  stdElemType elemType;
  ot::DA::index idx[8];

  int ***Aijk = (int ***)m_stencil;

  alignElementAndVertices(m_octDA, elemType, idx);       

  PetscScalar *fibers  = (PetscScalar *) m_farray;

  for (int k = 0;k < 8;k++) {
    for (int j=0;j<8;j++) {
      in[m_uiDof*idx[k]]   += fac*(Aijk[elemType][k][j]) * fibers[m_uiDof*idx[j]];   
      in[m_uiDof*idx[k]+1] += fac*(Aijk[elemType][k][j]) * fibers[m_uiDof*idx[j]+1]; 
      in[m_uiDof*idx[k]+2] += fac*(Aijk[elemType][k][j]) * fibers[m_uiDof*idx[j]+2]; 
    }//end for j
  }//end for k

#ifdef __DEBUG__  
  std::cout << "Leaving " << __func__ << std::endl;
#endif  
  return true;
}

#endif
