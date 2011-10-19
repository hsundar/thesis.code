#include "TransMatrix.h"
#include <float.h>


TransMatrix::TransMatrix()
{

}

TransMatrix::~TransMatrix()
{

}

TransMatrix::TransMatrix(double m01, double m02, double m03, double m04, 
		        double m05, double m06, double m07, double m08,
				double m09, double m10, double m11, double m12, 
				double m13, double m14, double m15, double m16 ) 
{
	m_MatElements[0] = m01; m_MatElements[1] = m02; m_MatElements[2] = m03; m_MatElements[3] = m04;  
	m_MatElements[4] = m05; m_MatElements[5] = m06; m_MatElements[6] = m07; m_MatElements[7] = m08;  
	m_MatElements[8] = m09; m_MatElements[9] = m10; m_MatElements[10] = m11; m_MatElements[11] = m12;  
	m_MatElements[12] = m13; m_MatElements[13] = m14; m_MatElements[14] = m15; m_MatElements[15] = m16;  
}

TransMatrix TransMatrix::inverse()
{
	TransMatrix matOut;
	matOut.SetElement(15,1.0);
	for (int i = 0; i < 3; i++) {
		matOut.SetElement(12+i, 0);
		for (int j = 0; j < 3; j++) {
			matOut.SetElement(4*i+j, GetElement(4*j+i));
			matOut.SetElement(12+i, matOut.GetElement(12+i) - GetElement(4*i+j) * GetElement(12+j));
		}
	}
	for (int i = 0; i < 4; i++) 
		matOut.SetElement(3+4*i,  GetElement(3+4*i));

	return matOut;
}

void TransMatrix::InvertMatrix(double* matOut)
{
	matOut[15]=1;
	for (int i = 0; i < 3; i++) {
		matOut[12+i]=0;
		for (int j = 0; j < 3; j++) {
			matOut[4*i+j]= GetElement(4*j+i);
			matOut[12+i]=  matOut[12+i] - GetElement(4*i+j) * GetElement(12+j);
		}
	}
}
void TransMatrix::Init(double *transArray)
{
	for (int i=0; i<16; i++){
		m_MatElements[i] = transArray[i];
	}
}

TransMatrix TransMatrix::MultMatrixLeftBy(TransMatrix& matIn)
{
	TransMatrix matOut;

	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += matIn.GetElement(k*4+j) * GetElement(i*4+k);
		}
	}
	for (int i = 0; i < 16; i++) 
		matOut.SetElement(i, tmp[i]);

	return matOut;
}
TransMatrix TransMatrix::MultMatrixLeftBy(double* matIn)
{
	TransMatrix matOut;

	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += matIn[k*4+j] * GetElement(i*4+k);
		}
	}
	for (int i = 0; i < 16; i++) 
		matOut.SetElement(i, tmp[i]);

	return matOut;
}
void TransMatrix::MultMatrixLeftBy(double* matIn, double* matOut)
{
	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += matIn[k*4+j] * GetElement(i*4+k);
		}
	}

	for (int i = 0; i < 16; i++) 
		matOut[i]=tmp[i];
}
TransMatrix TransMatrix::MultMatrixRightBy(TransMatrix& matIn)
{
	TransMatrix matOut;

	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += GetElement(k*4+j) * matIn.GetElement(i*4+k);
		}
	}
	for (int i = 0; i < 16; i++) 
		matOut.SetElement(i, tmp[i]);

	return matOut;
}

TransMatrix TransMatrix::MultMatrixRightBy(double* matIn)
{
	TransMatrix matOut;

	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += GetElement(k*4+j) * matIn[i*4+k];
		}
	}
	for (int i = 0; i < 16; i++) 
		matOut.SetElement(i, tmp[i]);

	return matOut;
}
void TransMatrix::MultMatrixRightBy(double* matIn, double* matOut)
{
	double tmp[16];

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++) {
			tmp[i*4+j] = 0;
			for (int k = 0; k < 4; k++) 
				tmp[i*4+j] += GetElement(k*4+j) * matIn[i*4+k];
		}
	}

	for (int i = 0; i < 16; i++) 
		matOut[i]=tmp[i];
}
void TransMatrix::CopyElementsTo(double* matOut)
{
	for (int i = 0; i < 16; i++) 
		matOut[i]=m_MatElements[i];
}

void TransMatrix::InvertMatrixWithScale(double* matOut)
{
	double det = m_MatElements[0]*(m_MatElements[5]*m_MatElements[10]-m_MatElements[9]*m_MatElements[6])-
				 m_MatElements[4]*(m_MatElements[1]*m_MatElements[10]-m_MatElements[9]*m_MatElements[2])+ 
				 m_MatElements[8]*(m_MatElements[1]*m_MatElements[6]-m_MatElements[5]*m_MatElements[2]);


	det = det + FLT_EPSILON;

	
	matOut[0] = ( m_MatElements[10] * m_MatElements[5] - m_MatElements[9] * m_MatElements[6])/det;
	matOut[1] = ( m_MatElements[9] * m_MatElements[2] - m_MatElements[10] * m_MatElements[1])/det;
	matOut[2] = ( m_MatElements[1] * m_MatElements[6] - m_MatElements[5] * m_MatElements[2])/det;
	matOut[3] = 0 ;

	matOut[4] = ( m_MatElements[8] * m_MatElements[6] - m_MatElements[10] * m_MatElements[4])/det;
	matOut[5] = ( m_MatElements[10] * m_MatElements[0] - m_MatElements[8] * m_MatElements[2])/det;
	matOut[6] = ( m_MatElements[4] * m_MatElements[2] - m_MatElements[0] * m_MatElements[6])/det;
	matOut[7] = 0;

	matOut[8] = ( m_MatElements[4] * m_MatElements[9] - m_MatElements[8] * m_MatElements[5])/det;
	matOut[9] = ( m_MatElements[8] * m_MatElements[1] - m_MatElements[0] * m_MatElements[9])/det;
	matOut[10]= ( m_MatElements[0] * m_MatElements[5] - m_MatElements[4] * m_MatElements[1])/det;
	matOut[11]= 0;

	matOut[12] = -(matOut[0] * m_MatElements[12] + matOut[4] * m_MatElements[13] + matOut[8 ] * m_MatElements[14] );
	matOut[13] = -(matOut[1] * m_MatElements[12] + matOut[5] * m_MatElements[13] + matOut[9 ] * m_MatElements[14] );
	matOut[14] = -(matOut[2] * m_MatElements[12] + matOut[6] * m_MatElements[13] + matOut[10] * m_MatElements[14] );
	matOut[15] = 1.0;
}

Point
TransMatrix::operator * (Point& ptIn) const
{
    Point ptOut;
	ptOut.x() = m_MatElements[0] * ptIn.x() + m_MatElements[4] * ptIn.y() + m_MatElements[ 8] * ptIn.z() + m_MatElements[12];
	ptOut.y() = m_MatElements[1] * ptIn.x() + m_MatElements[5] * ptIn.y() + m_MatElements[ 9] * ptIn.z() + m_MatElements[13];
	ptOut.z() = m_MatElements[2] * ptIn.x() + m_MatElements[6] * ptIn.y() + m_MatElements[10] * ptIn.z() + m_MatElements[14];

	return ptOut;
}
