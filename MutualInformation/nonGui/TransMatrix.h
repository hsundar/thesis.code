#ifndef _TRANSMATRIX_H
#define _TRANSMATRIX_H

#include "Point.h"

class TransMatrix  
{
public:

	/// Constructor
	TransMatrix();
	TransMatrix(double m01, double m02, double m03, double m04, 
		        double m05, double m06, double m07, double m08,
				double m09, double m10, double m11, double m12, 
		        double m13, double m14, double m15, double m16 );

	/// Destructor
	virtual ~TransMatrix();

	// static functions ...

	static TransMatrix Identity() { return TransMatrix(1.,0.,0.,0.,  0.,1.,0.,0.,  0.,0.,1.,0.,  0.,0.,0.,1.); } 
	static TransMatrix Translation(double x, double y, double z) { return TransMatrix(1.,0.,0.,0.,  0.,1.,0.,0.,  0.,0.,1.,0.,  x,y,z,1.); }
        static TransMatrix RotationX(double a) { return TransMatrix (1.,0.,0.,0., 0.,cos(a),-sin(a),0., 0.,sin(a),cos(a),0., 0.,0.,0.,1.);  }

	/**
	*  @brief		Initializes the Matrix by using a input matrix in a format of double[16]
	*  @param		transArray	:the input matrix in a format of double[16]
	**/
	void Init(double* transArray);

	/**
	*  @brief		gets this matrix elements in a array of double
	*  @return		the matrix elements in a format of double[16]
	**/
	double* GetDataPointer(){return m_MatElements;}

	/**
	*  @brief		Returns the element at indice i
	*  @param		i	: the element's indice in the matrix
	*  @return		the element in a type double
	**/
	double GetElement(int i){return m_MatElements[i];}

	/**
	*  @brief		Updates the matrix by writing a new value at indice i
	*  @param		i		:the new element's indice in the matrix
	*  @param		value	:the new value 
	**/
	void   SetElement(int i, double value){m_MatElements[i]= value;}

	/**
	*  @brief		Invert this matrix and returns the result in double[16] format
	*  @param		matOut	:the output matrix in double[16] format
	**/
	void InvertMatrix(double* matOut);

	/**
	*  @brief		Invert this matrix 
	*  @return		the inverted matrix in TransMatrix format
	**/
	TransMatrix inverse();

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from left (matIn x thisMatrix)  
	*				and returns the result in a TransMatrix format
	*  @param		matIn	:the input matrix in TransMatrix format
	*  @return		the resulting matrix in TransMatrix format
	**/
	TransMatrix MultMatrixLeftBy(TransMatrix& matIn);

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from left (matIn x thisMatrix)  
	*				and returns the result in a TransMatrix format
	*  @param		matIn	:the input matrix in a format of double[16]
	*  @return		the resulting matrix in TransMatrix format
	**/
	TransMatrix MultMatrixLeftBy(double* matIn);

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from left (matIn x thisMatrix)  
	*				and returns the result in matOut in double[16] format
	*  @param		matIn	:the input matrix in double[16] format
	*  @param		matOut  :the resulting matrix in double[16] format
	**/
	void MultMatrixLeftBy(double* matIn, double *matOut);

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from right (thisMatrix x matIn ) and returns 
	*				the result in a TransMatrix format
	*  @param		matIn	:the input matrix in TransMatrix format
	*  @return		the resulting matrix in TransMatrix format
	**/
	TransMatrix MultMatrixRightBy(TransMatrix& matIn);

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from right (thisMatrix x matIn) and returns 
	*				the result in matOut in double[16] format
	*  @param		matIn	:the input matrix in a format of double[16]
	*  @param		matOut  :the resulting matrix in double[16] format
	**/
	void MultMatrixRightBy(double* matIn, double *matOut);

	/**
	*  @brief		Multiplies this matrix by the input matrix matIn from right (thisMatrix x matIn) and returns 
	*				the result in a TransMatrix format
	*  @param		matIn	:the input matrix in double[16] format
	*  @return		the resulting matrix in TransMatrix format
	**/
	TransMatrix MultMatrixRightBy(double* matIn);

	/**
	*  @brief		copies elements of this matrix into matOut
	*  @param		matOut  :the outPut matrix in double[16] format
	**/
	void CopyElementsTo(double* matOut);

	/**
	*  @brief		inverts the matrix 
	*  @param		matOut  :the outPut matrix in double[16] format
	**/
	void InvertMatrixWithScale(double* matOut);

	/// Apply the transformation matrix to the given input point
	Point operator * ( Point& ptIn) const;

private:

	double m_MatElements[16]; ///< the matrix elements, unwrapped column-wise

};

#endif 
