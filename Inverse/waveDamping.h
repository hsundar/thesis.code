/**
 *  @file	waveDamping.h
 *  @brief	Main class to apply wave damping.
 *  @author Hari Sundar
 *  @date	  5/15/7
 * 
 *  Main class to apply wave damping.
 */


#ifndef __WAVE_DAMPING_H_
#define __WAVE_DAMPING_H_

#include "feMatrix.h"
#include "stiffnessMatrix.h"
#include "massMatrix.h"

/**
 *  @file	waveDamping.h
 *  @brief	Main class to apply wave damping.
 *  @author Hari Sundar
 *  @date	  5/15/7
 * 
 *  Main class to apply wave damping.
 */

class waveDamping : public feMatrix<waveDamping>
{
  public:	
    waveDamping(daType da);
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

    inline bool ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale);
    inline bool ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale);

    bool preMatVec();
    bool postMatVec();

    inline bool initStencils() {
      return true;
    }

    void setAlpha(double a) { 
      m_dAlpha = a;
    }

    void setBeta(double b) {
      m_dBeta = b;
    }

    void setMassMatrix (massMatrix * mass) {
      m_matMass = mass;
    }

    void setStiffnessMatrix (stiffnessMatrix *stiff) {
      m_matStiffness = stiff;
    }

   private:
    double              m_dAlpha;
    double              m_dBeta;

    massMatrix *          m_matMass;
    stiffnessMatrix *     m_matStiffness;
};

waveDamping::waveDamping(daType da) {
#ifdef __DEBUG__
  assert ( ( da == PETSC ) || ( da == OCT ) );
#endif
  m_daType = da;
  m_DA    = NULL;
  m_octDA   = NULL;
  m_stencil = NULL;
}

bool waveDamping::preMatVec() {
  m_matMass->preMatVec(); 
  m_matStiffness->preMatVec();
  // std::cout << "Leaving Damping preMatVec" << std::endl;
  return true;
}

bool waveDamping::postMatVec() {
  m_matMass->postMatVec();
  m_matStiffness->postMatVec();
  // std::cout << "Leaving Damping postMatVec" << std::endl;
  return true;
}

bool waveDamping::ElementalMatVec(int i, int j, int k, PetscScalar ***in, PetscScalar ***out, double scale) {
  return ( m_matMass->ElementalMatVec(i,j,k,in,out,scale*m_dAlpha) + 
           m_matStiffness->ElementalMatVec(i,j,k,in,out,scale*m_dBeta) );
}

bool waveDamping::ElementalMatVec(unsigned int idx, PetscScalar *in, PetscScalar *out, double scale) {
  return ( m_matMass->ElementalMatVec(idx, in, out, scale*m_dAlpha) +
           m_matStiffness->ElementalMatVec(idx,in,out,scale*m_dBeta) );
}

bool waveDamping::ElementalMatGetDiagonal(int i, int j, int k, PetscScalar ***diag, double scale) {
  return ( m_matMass->ElementalMatGetDiagonal(i,j,k,diag,scale*m_dAlpha) + 
           m_matStiffness->ElementalMatGetDiagonal(i,j,k,diag,scale*m_dBeta) );
}

bool waveDamping::ElementalMatGetDiagonal(unsigned int idx, PetscScalar *diag, double scale) {
  return ( m_matMass->ElementalMatGetDiagonal(idx, diag, scale*m_dAlpha) +
           m_matStiffness->ElementalMatGetDiagonal(idx,diag,scale*m_dBeta) );
}

#endif
