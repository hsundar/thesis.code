/**
 * @file Fields.h
 * @brief Contains 3 container classes (scalar, vector, tensor) derived from Field
*/

#include "Field.h"
//#include "Ellipsoid.h" 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <math.h>

/**
 * @brief a data structure to store a 3-d point
 */
struct Pt3d{
  float x,y,z;
};

/**
 * @brief A scalar field class 
*/
class ScalarField : public Field<float, 1>
{
 public:
  /** 
   * computes truncated Gaussian smoothing filter from input sigma parameter
   * filter itself is represented as a scalar field
   */
  void ComputeSmoothingFilter(float);

  void ComputeSmoothingFilter(float sigma, float cen_x, float cen_y, float cen_z);

  /** 
   * smooths a scalar field using the specified scalar field filter
   */
  void Smooth(ScalarField *smoothingFilter, ScalarField *out);
};

/**
 * @brief  A 3-d vector field class 
*/
class VectorField : public Field<float, 3>
{
 public:
  /**
   * Adds (x,y,z) offset to the vector field 
   */
  void AddOffset();

  /**
   * Subtracts (x,y,z) offset from the vector field 
   */
  void SubtractOffset();

  /**
   * Flips particular component of vector field 
   */
  void FlipComp(int);

  /**
   * Normalizes vector field to have unit length at all voxels
   */
  void Normalize();

  /**
   * Returns vector at particular position in vector field as a 3d point
   */
  int getVectorAtAnyPosition(float px, float py, float pz, struct Pt3d & P);

  /**
   * Reverses vector field
   */
  void reverseWarpField(VectorField *out);

   /**
   * Reverses vector field
   */
  void reverseWarpFieldHAMMER(VectorField *out);

  /** 
   * smooths vector field using the specified scalar field filter
   */
  void Smooth(ScalarField *smoothingFilter, VectorField *out);
};

class TensorField;


/**
 *@brief  A tensor field class of symmetric matrices represented as 6-d vectors
 * note that upon occasion, symmetric but not positive-definite matrices
 *  may also be stored in objects of the TensorField class 
*/
class TensorField : public Field<float, 6>
{
 public:
  /**
   * creates tensor field from 6 input scalar fields
   */
  void createInterleavedDTI(ScalarField *Dxx, ScalarField *Dyy, ScalarField *Dzz, ScalarField *Dxy, ScalarField *Dxz, ScalarField *Dyz);

  /**
   * extracts 6 scalar fields from tensor field
   */
  void extractDTIcomponents(ScalarField *Dxx, ScalarField *Dyy, ScalarField *Dzz, ScalarField *Dxy, ScalarField *Dxz, ScalarField *Dyz);

  /**
   * computes scalar Trace field from tensor field
   */
  void ComputeTrace(ScalarField *traceField);

  /**
   * computes scalar FA (fractional anisotropy) field from tensor field
   */
  void ComputeFA(ScalarField *faVol);

  /**
   * computes vector field corresponding to specified eigen vector of tensor field
   * options to specify which eigen-vector and whether vector field to weight by FA
   */
  void ComputePD(VectorField *PDVol, int d_index, int weight_flag);

  /**
   * computes 3 eigenvalue, 3 eigenvector fields from tensor field
   */
  void ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field);

  /**
   * computes 3 eigenvalue, 3 eigenvector, trace, FA fields from tensor field
   * there is a option to weight eigenvector fields by FA
   */
  void ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field, ScalarField *traceField,ScalarField *faField, int FAweighting);

  /**
   * computes 3 eigenvalue, 3 eigenvector, trace, FA, 3 Euler angle fields from tensor field
   * there is a option to weight eigenvector fields by FA
   */
  void ComputeEigD(ScalarField *e1Field,ScalarField *e2Field,ScalarField *e3Field,VectorField *PD1Field,VectorField *PD2Field,VectorField *PD3Field, ScalarField *traceField,ScalarField *faField, int FAweighting, ScalarField *PhiField, ScalarField *ThetaField, ScalarField *PsiField);

  /**
   * converts Euler angles to a rotation matrix 
   */
  void EulerAngles2RotMat(double *RotMat, float phi, float theta, float psi);

  /**
   * converts rotation matrix to Euler angles (refer Greg Slabaugh's tutorial)
   */
  void RotMat2EulerAngles(const double *RotMat, float &phi, float &theta, float &psi);

  /** 
   * smooths tensor field using the specified scalar field filter in Log-Euclidean domain
   */
  void logESmooth(ScalarField *smoothingFilter, TensorField *dti_out);

  /** 
   * smooths tensor field using the specified scalar field filter
   */
  void ESmooth(ScalarField *smoothingFilter, TensorField *dti_out);

  /** 
   * smooths tensor field using the specified scalar field filter within a mask  of non-zero tensors
   */
  void ESmooth(ScalarField *smoothingFilter, TensorField *dti_out, unsigned char *mask) ;

  /**
   * computes symmetric matrix logarithm of tensor field within a mask of non-zero tensors
   */
  void logTField(TensorField *dti_out, unsigned char *mask);

  /**
   * computes matrix exponential of input symmetric matrix field
   */
  void expTField(TensorField *dti_out, unsigned char *mask);

  /**
   * computes non-zero tensor mask using a trace-based threshold
   * a threshold factor is applied to the 
   * the actual trace-based threshold is saved in the public variable traceTensorThreshold
   */
  void computeNZmask(unsigned char *mask, float threshold_factor); 

  float traceTensorThreshold;  // a public variable to store trace threshold fordeciding non-zero tensor

  /**
   * a function to re-orient the tensor field using input deformation field
   * the finite-strain approximation to compute re-orientation matrix
   */
  void reOrientTensorFieldFS(VectorField *warpField, TensorField *dti_out);

  // these are two slow, but accurate functions
  /**
   * displace the tensor field using input forward deformation field and sigma specifying kernel for Gaussian interpolation
   */
  void displaceTensorFieldUsingFwdField(VectorField *warpField, TensorField *dti_out, float sigma, float filter_size_red_factor);

  /**
   * displace the tensor field using input inverse deformation field and sigma specifying kernel for Gaussian interpolation
   */
  void displaceTensorFieldUsingInvField(VectorField *warpField, TensorField *dti_out, float sigma, float filter_size_red_factor);

  // this is a faster function
  /**
   * displace the tensor field using input inverse deformation field and trilinear interpolation
   */
  void displaceTensorFieldUsingInvFieldTL(VectorField *warpField, TensorField *dti_out);
  
 
};


