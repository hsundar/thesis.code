/**
*  @file	SimilarityMeasures.h
*  @brief	Computes various similarity measures for 2-dimensional images.
*  @author	Hari Sundar
*  @date	9/13/2003
*
*
*  So far those measures are only tested using the template type unsigned char.
*  Using more precise data may result in overflow of certain sums, as I am still
*  using long variables for summation due to higher speed. Changing all long
*  variables to double should solve this problem. Note that the gradient images
*  are saved as int, i.e. 32 bit signed values. Therefore it makes no sense to
*  use gradient-based measures with template data of a higher precision.
*
*  Stereo Images are supported, i.e. two images appended horizontally. Using the
*  respective mode with SetStereo(true) takes care of excluding the border in the
*  middle of the combined image out of the measure computation.
*
*  A Region of Interest may be specified with SetMask(mask) and SetUseMask(true),
*  mask is unsigned char ftk image, with all non-zero values regarded as included
*  in the ROI.
*
*  Before the class instance is initialized with InitInternalBuffer(...), both
*  images have to be supplied with SetFirstImage(...) and SetSecondImage(...) in
*  order to have the size information. The image data may be changed later, though.
**/


#ifndef SBIA_SIMILARITYMEASURES_H
#define SBIA_SIMILARITYMEASURES_H

#include <vector>

#include "array1d.h"
#include "array2d.h"
#include "image/Image.h"
#include "image/IteratorPosition.h"
#include "image/Iterator.h"

template <typename T>
class SimilarityMeasures
{
public:

	/**
	*  @brief		Constructor
	*
	*  Resets all variables to their initial values, no allocation is done yet.
	*
	**/
	SimilarityMeasures();


	/**
	*  @brief		Destructor
	*
	*  Frees all the allocated memory.
	*
	**/
	virtual ~SimilarityMeasures();


	/**
	*  @brief		Enables or disables the stereo mode
	*  @param		flag	true=stereo, false=mono
	*
	*  In stereo mode, the images consist of two parts, appended
	*  horizontally. Therefore the gradient- and neighborhood-based
	*  measures have to omit some columns in the middle of the image,
	*  where the border between the two single images lies.
	*
	**/
	void SetStereo(bool flag) {m_stereo = flag;};


	/**
	*  @brief		Enables or disables the ROI mask.
	*  @param		flag	true=on, false=off
	*
	*  If enabled, a binary region of interest mask has to be supplied
	*  earlier with the SetMask method.
	*
	**/
	void SetUseMask(bool flag) {m_useMask = flag;};


	/**
	*  @brief		Set the intensity threshold
	*  @param		value	the threshold value
	*
	*  If set to any value greater zero, all pixels whose intensity is smaller
	*  are being ignored in the measure computation. This is so far only implemented
	*  for the histogram-based measures, MI and CR. Further extension is probably
	*  not necessary, as this is only for testing purposes and any region of interest
	*  may be defined via the binary mask anyway.
	*
	**/
	void SetIgnoreIntensitiesBelow(T value) {m_ignoreIntensitiesBelow = value;}


	/**
	*  @brief		To initialize the internal buffers.
	*  @param		Number of bits for the histogram.
	*  @return		False if there is a problem allocating the buffers.
	**/
	bool InitInternalBuffer(int histBits = 8);


	/**
	*  @brief		Computes the measures.
	*  @return		returns false if some error happens.
	*
	*  The input has to be provided already, also the type of the measure should
	*  be specified in advance.
	**/
	bool CalculateMeasures();


	/** @name Specify the input data **/
	//@{
		void SetFirstImage(image::Image<T, 2>& rImage1, T image1MinValue, int image1BitsUsed);
		void SetSecondImage(image::Image<T, 2>& rImage2, T image2MinValue, int image2BitsUsed);
		void SetMask(image::Image<unsigned char, 2>& rMask);
	//@}


	/** @name Set the type of measures to be computed **/
	//@{
		void SetMutualInformation(bool flag) {m_calculateMI = flag;};
		void SetMutualInformationNorm(bool flag){m_normalizeMI = flag;};
		void SetSumofSquareDifference(bool flag){m_calculateSSD = flag;};
		void SetSumofAbsoluteDifference(bool flag){m_calculateSAD = flag;};
		void SetNormalizedCrossCorr(bool flag){m_calculateNCC = flag;};
		void SetLocalNormalizedCorr(bool flag){m_calculateLNC = flag;};
		void SetVarianceWeightedCorr(bool flag){m_calculateVWC = flag;};
		void SetPatternIntensity(bool flag){m_calculatePI = flag;};
		void SetGradientCorrelation(bool flag){m_calculateGC = flag;};
		void SetGradientDifference(bool flag){m_calculateGD = flag;};
		void SetCorrelationRatio(bool flag){m_calculateCR = flag;};
	//@}
	

	/** @name Access the results of the computation **/
	//@{
		double GetMutualInformation(){return m_similarityMI;};
		double GetSumofSquareDifference(){return m_similaritySSD;};
		double GetSumofAbsoluteDifference(){return m_similaritySAD;};
		double GetNormalizedCrossCorr(){return m_similarityNCC;};
		double GetLocalNormalizedCorr(){return m_similarityLNC;};
		double GetVarianceWeightedCorr(){return m_similarityVWC;};
		double GetPatternIntensity(){return m_similarityPI;};
		double GetGradientCorrelation(){return m_similarityGC;};
		double GetGradientDifference(){return m_similarityGD;};
		double GetCorrelationRatio(){return m_similarityCR;};
		double GetFirstEntropy(){return m_entropy1;};
		double GetSecondEntropy(){return m_entropy2;};
		double GetJointEntropy(){return m_entropy3;};
	//@}


	/** @name Access histogram information */
	//@{
		hs::array1d<unsigned int>* GetFirstHistogram() { return m_phist1; }
		hs::array1d<unsigned int>* GetSecondHistogram() { return m_phist2; }
		hs::array2d<unsigned int>* GetJointHistogram() { return m_pDistribution; }
	//@}


	double m_lncError;									 // for validation of LNC: occurrences of zero variance

private:

	// Cut out of the main routine: computes Gradient Correlation for existing gradient images
	double CalculateGradientCorrelation();

	// The same for Gradient Difference
	double CalculateGradientDifference();

	/** @name Results of the measure computation */
	//@{
		double m_similarityMI;		// Mutual Information
		double m_similaritySSD;		// Sum of Squared Differences
		double m_similaritySAD;		// Sum of Absolute Differences
		double m_similarityNCC;		// Normalized Cross Correlation
		double m_similarityLNC;		// Local Normalized Correlation
		double m_similarityVWC;		// Variance Weighted Local Normalized Correlation
		double m_similarityPI;      // Pattern Intensity
		double m_similarityGC;		// Gradient Correlation
		double m_similarityGD;		// Gradient Difference
		double m_similarityCR;		// Correlation Ratio
		double m_entropy1;			// Entropy of the first image
		double m_entropy2;			// Entropy of the second image
		double m_entropy3;			// Joint entropy of both images
	//@}

	/** @name Flags which measures to compute */
	//@{
		bool  m_calculateMI;
		bool  m_calculateSSD;
		bool  m_calculateSAD;
		bool  m_calculateNCC;
		bool  m_calculateLNC;
		bool  m_calculateVWC;
		bool  m_calculatePI;
		bool  m_calculateGC;
		bool  m_calculateGD;
		bool  m_calculateCR;
	//@}

	/** @name Variances of the gradient images */
	//@{
		double m_varH1;
		double m_varH2;
		double m_varV1;
		double m_varV2;
	//@}

	bool m_normalizeMI;									 // compute Mutual Information normalized or not
	bool m_stereo;										 // use two images attached to each other
	bool m_useMask;										 // we are using a binary ROI mask
	T m_image1MinValue;									 // lowest value in first image
	T m_image2MinValue;									 // lowest value in second image
	int m_image1BitsUsed;								 // number of used bits in first image
	int m_image2BitsUsed;								 // number of used bits in second image
	hs::array1d<unsigned int>*  m_phist1;				 // pointer to the first histogram
	hs::array1d<unsigned int>*  m_phist2;				 // pointer to the second histogram
	hs::array2d<unsigned int>*  m_pDistribution;		 // 2D integer array containing the distribution
	hs::image::Image<T, 2>*  m_pImage1;				 // data of first (moving) image
	hs::image::Image<T, 2>*  m_pImage2;				 // data of second (fixed) image
	hs::image::Image<int, 2>*		     m_pGradientH1;	 // horizontal gradient of first image
	hs::image::Image<int, 2>*		     m_pGradientH2;	 // horizontal gradient of second image
	hs::image::Image<int, 2>*           m_pGradientV1;	 // vertical gradient of first image
	hs::image::Image<int, 2>*           m_pGradientV2;	 // vertical gradient of second image
	hs::image::Image<int, 2>*           m_pDiffImage;	 // difference image for PI
	hs::image::Image<unsigned char,2>*	 m_pMask;		 // binary ROI mask
	hs::image::Image<unsigned char,2>*	 m_pMaskGradient;// eroded mask for gradient calculation
	int m_histBits;										 // number of bits used for histograms
	double m_sigmaPI;   								 // square of sigma for pattern intensity
	int m_radiusPI;										 // radius for pattern intensity
	int m_windowLNC;									 // window size for LNC (half of the width!)
	int m_gradSize;										 // number of used pixels in the gradient images
	T m_ignoreIntensitiesBelow;							 // for testing strategies to omit dark image regions
};

/* ----------- End of class definition. Here comes the implementation. ---------- */


template <typename T>
SimilarityMeasures<T>::SimilarityMeasures()
{
	m_phist1 = NULL;
	m_phist2 = NULL;
	m_pDistribution = NULL;
	m_pGradientH1 = NULL;
	m_pGradientH2 = NULL;
	m_pGradientV1 = NULL;
	m_pGradientV2 = NULL;
	m_pDiffImage  = NULL;
	m_pMask       = NULL;
	m_pMaskGradient = NULL;
	m_histBits = 0;
	m_sigmaPI  = 100;
	m_radiusPI = 3;
	m_windowLNC = 3;
	m_gradSize = 0;
	m_lncError = 0;
	m_ignoreIntensitiesBelow = 0;
	m_similarityMI  = 0.0;
	m_similaritySSD = 0.0;	
	m_similaritySAD = 0.0;	
	m_similarityNCC = 0.0;
	m_similarityLNC = 0.0;
	m_similarityVWC = 0.0;
	m_similarityPI  = 0.0; 
	m_similarityGC  = 0.0;
	m_similarityGD  = 0.0;
	m_similarityCR  = 0.0;
	m_entropy1      = 0.0;
	m_entropy2      = 0.0;
	m_entropy3      = 0.0;
	m_calculateMI  = true;
	m_calculateSSD = false;
	m_calculateSAD = false;
	m_calculateNCC = false;
	m_calculateLNC = false;
	m_calculateVWC = false;
	m_calculatePI  = false;
	m_calculateGC  = false;
	m_calculateGD  = false;
	m_calculateCR  = false;
	m_normalizeMI  = true;
	m_stereo       = false;
	m_useMask      = false;
}


template <typename T>
SimilarityMeasures<T>::~SimilarityMeasures()
{
	if (m_phist1) delete m_phist1;
	if (m_phist2) delete m_phist2;
	if (m_pDistribution) delete m_pDistribution;
	if (m_pGradientH1) delete m_pGradientH1;
	if (m_pGradientH2) delete m_pGradientH2;
	if (m_pGradientV1) delete m_pGradientV1;
	if (m_pGradientV2) delete m_pGradientV2;
	if (m_pDiffImage) delete m_pDiffImage;
	if (m_pMaskGradient) delete m_pMaskGradient;
}


template <typename T>
void SimilarityMeasures<T>::SetFirstImage(image::Image<T, 2>& rImage1, T image1MinValue, int image1BitsUsed)
{
	m_pImage1 = &rImage1;
	m_image1MinValue =  image1MinValue;
	m_image1BitsUsed =  image1BitsUsed;
}


template <typename T>
void SimilarityMeasures<T>::SetSecondImage(image::Image<T, 2>& rImage2, T image2MinValue, int image2BitsUsed)
{
	m_pImage2 = &rImage2;
	m_image2MinValue =  image2MinValue;
	m_image2BitsUsed =  image2BitsUsed;
}


template <typename T>
void SimilarityMeasures<T>::SetMask(image::Image<unsigned char, 2>& rMask) {
	m_pMask = &rMask;
	if (m_pMaskGradient) delete m_pMaskGradient;
	hs::Size<hs::IndexType, 2> size( m_pImage1->getSize(0)-2, m_pImage1->getSize(1)-2 );
	m_pMaskGradient = hs::image::Image<unsigned char, 2> (size);
	// Perform the erosion for the gradient mask
	for (int y = 1; y < m_pImage1->getSize(1)-1; y++) {
		for (int x = 1; x < m_pImage1->getSize(0)-1; x++) {
			(*m_pMaskGradient)(x-1, y-1) = (
				((*m_pMask)(x-1,y-1)) &&
				((*m_pMask)(x  ,y-1)) &&
				((*m_pMask)(x+1,y-1)) &&
				((*m_pMask)(x-1,y  )) &&
				((*m_pMask)(x  ,y  )) &&
				((*m_pMask)(x+1,y  )) &&
				((*m_pMask)(x-1,y+1)) &&
				((*m_pMask)(x  ,y+1)) &&
				((*m_pMask)(x+1,y+1))
			);
		}
	}
}


template <typename T>
bool SimilarityMeasures<T>::InitInternalBuffer(int histBits)
{
	// delete existing buffers
	if (m_phist1) { delete m_phist1; m_phist1 = NULL; }
	if (m_phist2) { delete m_phist2; m_phist2 = NULL; }
	if (m_pDistribution) { delete m_pDistribution; m_pDistribution = NULL; }
	if (m_pGradientH1) { delete m_pGradientH1; m_pGradientH1 = NULL; }
	if (m_pGradientV1) { delete m_pGradientV1; m_pGradientV1 = NULL; }
	if (m_pGradientH2) { delete m_pGradientH2; m_pGradientH2 = NULL; }
	if (m_pGradientV2) { delete m_pGradientV2; m_pGradientV2 = NULL; }
	if (m_pDiffImage) { delete m_pDiffImage; m_pDiffImage = NULL; }
	
	// buffers for histogram based measures
	if ((m_calculateMI) || (m_calculateCR)) {
		m_phist1 = new hs::array1d<unsigned int>(1 << histBits);
		if (m_phist1 == NULL) return false;

		m_phist2 = new hs::array1d<unsigned int>(1 << histBits);
		if (m_phist2 == NULL) return false;

		m_pDistribution = new hs::array2d<unsigned int>(1 << histBits,1 << histBits);
		if (m_pDistribution == NULL) return false;
	}

	// buffers for gradient based measures
	if ((m_calculateGC) || (m_calculateGD)) {
		hs::Size<hs::IndexType, 2> size1( m_pImage1->getSize(0)-2, m_pImage1->getSize(1)-2 );
		hs::Size<hs::IndexType, 2> size2( m_pImage2->getSize(0)-2, m_pImage2->getSize(1)-2 );
		m_pGradientH1 = new hs::image::Image<int, 2> (size1);
		if (m_pGradientH1 == NULL) return false;

		m_pGradientV1 = new hs::image::Image<int, 2> (size1);
		if (m_pGradientV1 == NULL) return false;

		m_pGradientH2 = new hs::image::Image<int, 2> (size2);
		if (m_pGradientH2 == NULL) return false;

		m_pGradientV2 = new hs::image::Image<int, 2> (size2);
		if (m_pGradientV2 == NULL) return false;
	}
	
	if (m_calculatePI) {
		hs::Size<hs::IndexType, 2> size( m_pImage1->getSize(0), m_pImage1->getSize(1));
		m_pDiffImage = new hs::image::Image<int, 2> (size);
	}

	m_histBits = histBits;

	return true;
}


template <typename T>
bool SimilarityMeasures<T>::CalculateMeasures()
{
	// reset the similarity measures
	m_similarityMI  = 0.0;
	m_similaritySSD = 0.0;
	m_similaritySAD = 0.0;
	m_similarityNCC = 0.0;
	m_similarityLNC = 0.0;
	m_similarityPI  = 0.0;
	m_similarityGC  = 0.0;
	m_similarityGD  = 0.0;
	m_similarityCR  = 0.0;
	m_entropy1      = 0.0;
	m_entropy2      = 0.0;
	m_entropy3      = 0.0;
	m_gradSize = 0;

	// clear histograms
	if ((m_calculateMI) || (m_calculateCR)) {
		m_phist1->fill(0);
		m_phist2->fill(0);
		m_pDistribution->fill(0);
	}

	int height = m_pImage1->getSize(1);
	int width  = m_pImage1->getSize(0);
	int size = 0;	// I compute the size by summing up the used pixels later

	// helper variables
	// I have no idea why long values are 32 bit here, same as int???
	// So I have to use double for some of the sums.
	long sum1 = 0, sum2 = 0, sum6 = 0, sum7 = 0;
	double sad = 0, ssd = 0, sum3 = 0, sum4 = 0, sum5 = 0;

	double countLNC = (2*m_windowLNC+1)*(2*m_windowLNC+1);
	double sumVWC = 0.0;

	// find the largest bit number used in the two images
	int bitsUsed = (m_image1BitsUsed > m_image2BitsUsed)? m_image1BitsUsed: m_image2BitsUsed;

	// first pass throught the images
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			T pixel1 = ((*m_pImage1)(x,y) -m_image1MinValue) << (bitsUsed - m_image1BitsUsed);
			T pixel2 = ((*m_pImage2)(x,y) -m_image2MinValue) << (bitsUsed - m_image2BitsUsed);
			int diff = pixel2 - pixel1;
			if (m_calculatePI) (*m_pDiffImage)(x,y) = diff;

			if ((!m_useMask) || ((*m_pMask)(x,y))) {
				sum1 += pixel1;
				sum2 += pixel2;
				if (m_calculateNCC) {
					sum3 += pixel1 * pixel1;
					sum4 += pixel2 * pixel2;
					sum5 += pixel1 * pixel2;
				}

				if (m_calculateSAD) sad += abs(diff);
				if (m_calculateSSD) ssd += diff * diff;

				// we need histograms for both MI and CR
				if ((m_calculateMI) || (m_calculateCR)) {
					// scale values for histogram indexing
					unsigned int histIndex1;
					unsigned int histIndex2;
					if (bitsUsed - m_histBits > 0) {
						histIndex1 = pixel1 >> (bitsUsed - m_histBits);
						histIndex2 = pixel2 >> (bitsUsed - m_histBits);
					}
					else {
						histIndex1 = pixel1 << (m_histBits - bitsUsed);
						histIndex2 = pixel2 << (m_histBits - bitsUsed);
					}

					if ((pixel1 >= m_ignoreIntensitiesBelow) &&
						(pixel2 >= m_ignoreIntensitiesBelow)) {
						// increase respective pixel in PDF
						m_pDistribution->fastGet(histIndex2, histIndex1)++;
						// increase histogram pixels
						m_phist1->fastGet(histIndex1)++;
						m_phist2->fastGet(histIndex2)++;
					}
				}
				size++;
			}

			// create gradient images
			if ((m_calculateGC) || (m_calculateGD)) {
				// we do not consider the border pixels
				if ((x) && (y) && (x < width-1) && (y < height-1)
					&& ((!m_stereo) || (x < width/2 - 1) || (x > width/2))
					&& ((!m_useMask) || ((*m_pMaskGradient)(x-1,y-1)))) {
					// we can neglect adding the minimum value, and shift the bits at the end
					(*m_pGradientH1)(x-1,y-1) = (
						-(*m_pImage1)(x-1,y-1)
						-2*(*m_pImage1)(x,y-1)
						-(*m_pImage1)(x+1,y-1)
						+(*m_pImage1)(x-1,y+1)
						+2*(*m_pImage1)(x,y+1)
						+(*m_pImage1)(x+1,y+1))
						<< (bitsUsed - m_image1BitsUsed);

					(*m_pGradientH2)(x-1,y-1) = (
						-(*m_pImage2)(x-1,y-1)
						-2*(*m_pImage2)(x,y-1)
						-(*m_pImage2)(x+1,y-1)
						+(*m_pImage2)(x-1,y+1)
						+2*(*m_pImage2)(x,y+1)
						+(*m_pImage2)(x+1,y+1))
						<< (bitsUsed - m_image2BitsUsed);

					(*m_pGradientV1)(x-1,y-1) = (
						-(*m_pImage1)(x-1,y-1)
						+(*m_pImage1)(x+1,y-1)
						-2*(*m_pImage1)(x-1,y)
						+2*(*m_pImage1)(x+1,y)
						-(*m_pImage1)(x-1,y+1)
						+(*m_pImage1)(x+1,y+1))
						<< (bitsUsed - m_image1BitsUsed);

					(*m_pGradientV2)(x-1,y-1) = (
						-(*m_pImage2)(x-1,y-1)
						+(*m_pImage2)(x+1,y-1)
						-2*(*m_pImage2)(x-1,y)
						+2*(*m_pImage2)(x+1,y)
						-(*m_pImage2)(x-1,y+1)
						+(*m_pImage2)(x+1,y+1))
						<< (bitsUsed - m_image2BitsUsed);
					m_gradSize++;
				}
			}

			// Local normalized correlation
			// Does not respect the bitmask yet, I still need to think about that...
			if ((m_calculateLNC) || (m_calculateVWC)) {
				if ((x >= m_windowLNC) && (x < width - m_windowLNC)
				&& (y >= m_windowLNC) && (y < height - m_windowLNC)
				&& ((!m_stereo) || (x < width/2 - m_windowLNC) || (x >= width/2 + m_windowLNC))) {
					int s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0;
					for (int yy = y - m_windowLNC; yy <= y + m_windowLNC; yy++) {
						for (int xx = x - m_windowLNC; xx <= x + m_windowLNC; xx++) {
							int p1 = ((*m_pImage1)(xx,yy) -m_image1MinValue) << (bitsUsed - m_image1BitsUsed);
							int p2 = ((*m_pImage2)(xx,yy) -m_image2MinValue) << (bitsUsed - m_image2BitsUsed);
							s1 += p1; s2 += p2;	s3 += p1 * p1; s4 += p2 * p2; s5 += p1 * p2;
						}
					}
					double co = ((double)s5) - ((double)(s1 * s2)) / countLNC;
					double v1 = ((double)s3) - ((double)(s1 * s1)) / countLNC;
					double v2 = ((double)s4) - ((double)(s2 * s2)) / countLNC;
					// We omit the window if one of the variances is zero, as the
					// correlation value is undefined then. See my thesis for a further
					// description why this is allowed :-)
					if ((v2 != 0) && (v1 != 0)) {
						if (m_calculateLNC) m_similarityLNC += co / (sqrt(v1) * sqrt(v2));
						if (m_calculateVWC) m_similarityVWC += v1 * co / (sqrt(v1) * sqrt(v2));
						sumVWC += v1;
						sum7++;
					}
				}
			}
		}
	}
	// end of first run through the images

	// Pattern Intensity
	// If we're talkin' about stereo images, we run over the
	// left and right part of them separately.
	if (m_calculatePI){
		for (int s = 0; s < (m_stereo? 2: 1); s++) {
			int w = m_stereo? width/2: width;
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < w; x++) {
					for (int yy = y - m_radiusPI; yy <= y + m_radiusPI; yy++) {
						for (int xx = x - m_radiusPI; xx <= x + m_radiusPI; xx++) {
							if ((xx >= 0) && (xx < w) && (yy >= 0) && (yy < height)
							&& ((x-xx)*(x-xx)+(y-yy)*(y-yy) <= m_radiusPI*m_radiusPI)) {
								if ((!m_useMask) || (((*m_pMask)(s*w+x,y)) && ((*m_pMask)(s*w+xx,yy)))) {
									double denom = (*m_pDiffImage)(s*w+x,y) - (*m_pDiffImage)(s*w+xx,yy);
									m_similarityPI += m_sigmaPI / (m_sigmaPI + denom*denom);
									sum6++;
								}
							}
						}
					}
				}
			}
		}
	}


	double mean1 = (double)sum1/(double)size;
	double mean2 = (double)sum2/(double)size;
	if (m_calculatePI)  m_similarityPI /= (double)sum6;
	if (m_calculateLNC) m_similarityLNC /= (double)sum7;
	if (m_calculateVWC) m_similarityVWC /= sumVWC;
	if (m_calculateSSD) m_similaritySSD = (double)ssd/(double)size;
	if (m_calculateSAD) m_similaritySAD = (double)sad/(double)size;

	if (m_calculateNCC) {
		double var1 = 0, var2 = 0, corr = 0;
/*		image::Iterator<T, 2> it1( *m_pImage1 );
		image::Iterator<T, 2> it2( *m_pImage2 );
		for (it1.goToBegin(); !it1.isAtEnd(); it1++, it2++) {
			T pixel1 = (*it1 -m_image1MinValue) << (bitsUsed - m_image1BitsUsed);
			T pixel2 = (*it2 -m_image2MinValue) << (bitsUsed - m_image2BitsUsed);
			corr += (pixel1-mean1)*(pixel2-mean2);
			var1 += (pixel1-mean1)*(pixel1-mean1);
			var2 += (pixel2-mean2)*(pixel2-mean2);
		} */

		corr = ((double)sum5) - ((double)(sum1 * sum2)) / (double)size;
		var1 = ((double)sum3) - ((double)(sum1 * sum1)) / (double)size;
		var2 = ((double)sum4) - ((double)(sum2 * sum2)) / (double)size;

		m_similarityNCC = corr / (sqrt(var1) * sqrt(var2));

	}

	int histSize = 1 << m_histBits;

	// calculate mutual information
	if (m_calculateMI){
		double p;
		for (int y = 0; y < histSize; y++) {
			p = ((double) m_phist1->fastGet(y)) / ((double) size);
			if (p > 0.0) m_entropy1 -= p * log(p);
			p = ((double) m_phist2->fastGet(y)) / ((double) size);
			if (p > 0.0) m_entropy2 -= p * log(p);
			for (int x = 0; x < histSize; x++) {
				p = ((double) m_pDistribution->fastGet(y,x)) / ((double) size);
				if (p > 0.0) m_entropy3 -= p * log(p);
			}
		}
		if (m_normalizeMI) {
			if (m_entropy1 + m_entropy2 > 0.0)
				m_similarityMI = 2 - 2 * m_entropy3 / (m_entropy1 + m_entropy2);
		}
		else
			m_similarityMI = m_entropy1 + m_entropy2 - m_entropy3;
	}

	/* calculate correlation ratio
	This measure is not symetric, therefore it makes a difference if the images
	are swapped! In this implementation, the first image is always the template image,
	i.e. it is seen as a base for the estimation of the second image. */
	if (m_calculateCR) {
		int i, j;
		// calculate variance and expectation out of the second histogram
		double var = 0, exp = 0;
		for (j = 0; j < histSize; j++) exp += j * m_phist2->fastGet(j);
		exp /= (double)size;
		for (j = 0; j < histSize; j++) var += j * j * m_phist2->fastGet(j);
		var = (var / (double)size) - exp * exp;
		// variance of expectation of conditional probability of a certain intensity :-)
		double v1 = 0, v2, px, ex;
		for (i = 0; i < histSize; i++) {
			ex = 0; v2 = 0;
			px = m_phist1->fastGet(i) / (double)size;
			if (px > 0.0) {
				for (j = 0; j < histSize; j++) ex += j * m_pDistribution->fastGet(j, i);
				ex /= (px * (double)size);
				for (j = 0; j < histSize; j++) v2 += j * j * m_pDistribution->fastGet(j, i);
				v2 = (v2 / (px * (double)size)) - ex * ex;
				
				v1 += v2 * m_phist1->fastGet(i);
			}
		}
		v1 /= (var * (double)size);
		m_similarityCR = 1 - v1;
	}

	// gradient difference needs values from the correlation, too
	if ((m_calculateGC) || (m_calculateGD)) m_similarityGC = CalculateGradientCorrelation();
	if (m_calculateGD) m_similarityGD = CalculateGradientDifference();


	return true;
	// This is a totally weird idea that I have to try. Rank correlation. Let's go...
	double R[256], S[256];
	double mean_r = 0, mean_s = 0, var1 = 0, var2 = 0, corr = 0;
	int sumh1 = 0, sumh2 = 0;
	// Omit black pixels
	for (int i = 1; i < histSize; i++) {
		int h1 = m_phist1->fastGet(i);
		int h2 = m_phist2->fastGet(i);
		R[i] = 0.5 * (double)(h1 + 1) + sumh1;
		S[i] = 0.5 * (double)(h2 + 1) + sumh2;
		sumh1 += h1; sumh2 += h2;
		mean_r += R[i] * (double)h1;
		mean_s += S[i] * (double)h2;
	}
	mean_r /= (double)histSize;
	mean_s /= (double)histSize;
	// And now run through the images!
	hs::image::Iterator<hs::image::Image<T, 2> > it1( *m_pImage1 );
	hs::image::Iterator<hs::image::Image<T, 2> > it2( *m_pImage2 );
	for (it1.goToBegin(); !it1.isAtEnd(); it1++, it2++) {
		T p1 = (*it1 - m_image1MinValue) << (bitsUsed - m_image1BitsUsed);
		T p2 = (*it2 - m_image2MinValue) << (bitsUsed - m_image2BitsUsed);
		if ((p1) && (p2)) {
			corr += (R[p1] - mean_r) * (S[p2] - mean_s);
			var1 += (R[p1] - mean_r) * (R[p1] - mean_r);
			var2 += (S[p2] - mean_s) * (S[p2] - mean_s);
		}
	}
	var1 = sqrt(var1) * sqrt(var2);
	m_similarityMI = 1.0 - (100000.0 * (var1 - corr) / var1);
	// End of this stupid experiment :-)

	return true;
}


template <typename T>
double SimilarityMeasures<T>::CalculateGradientCorrelation() {

	int width  = m_pGradientH1->getSize(0);
	int height = m_pGradientH1->getSize(1);
	// Same problem, long is 32 bit only.
	long sum1 = 0, sum2 = 0;
	double sum3 = 0, sum4 = 0, sum5 = 0;
	double corr, var1, var2;

	// set iterators
	image::IteratorPosition<image::Image<int, 2> > itH1( *m_pGradientH1 );
	image::IteratorPosition<image::Image<int, 2> > itH2( *m_pGradientH2 );
	image::IteratorPosition<image::Image<int, 2> > itV1( *m_pGradientV1 );
	image::IteratorPosition<image::Image<int, 2> > itV2( *m_pGradientV2 );

	// correlation of horizontal gradients
	for (; !itH1.isAtEnd(); itH1++, itH2++) {
		int x = itH1.getIndex()[0];
		int y = itH1.getIndex()[1];
		if (((!m_stereo) || (x < width/2 - 1) || (x > width/2))
		&& ((!m_useMask) || ((*m_pMaskGradient)(x,y)))) {
			sum1 += *itH1;
			sum2 += *itH2;
			sum3 += (*itH1) * (*itH1);
			sum4 += (*itH2) * (*itH2);
			sum5 += (*itH1) * (*itH2);
		}
	}
	corr = ((double)sum5) - ((double)(sum1 * sum2)) / (double)m_gradSize;
	var1 = ((double)sum3) - ((double)(sum1 * sum1)) / (double)m_gradSize;
	var2 = ((double)sum4) - ((double)(sum2 * sum2)) / (double)m_gradSize;
	m_varH1 = var1 / m_gradSize; m_varH2 = var2 / m_gradSize;
	double nccH = corr / (sqrt(var1) * sqrt(var2));
	sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0;

	// correlation of vertical gradients
	for (; !itV1.isAtEnd(); itV1++, itV2++) {
		int x = itV1.getIndex()[0];
		int y = itV1.getIndex()[1];
		if (((!m_stereo) || (x < width/2 - 1) || (x > width/2))
		&& ((!m_useMask) || ((*m_pMaskGradient)(x,y)))) {
			sum1 += *itV1;
			sum2 += *itV2;
			sum3 += (*itV1) * (*itV1);
			sum4 += (*itV2) * (*itV2);
			sum5 += (*itV1) * (*itV2);
		}
	}
	corr = ((double)sum5) - ((double)(sum1 * sum2)) / (double)m_gradSize;
	var1 = ((double)sum3) - ((double)(sum1 * sum1)) / (double)m_gradSize;
	var2 = ((double)sum4) - ((double)(sum2 * sum2)) / (double)m_gradSize;
	m_varV1 = var1 / m_gradSize; m_varV2 = var2 / m_gradSize;
	double nccV = corr / (sqrt(var1) * sqrt(var2));

	return (nccH + nccV) / 2;
}


template <typename T>
double SimilarityMeasures<T>::CalculateGradientDifference() {
	int width  = m_pGradientH1->getSize(0);
	int height = m_pGradientH1->getSize(1);
	double sum = 0;
	image::IteratorPosition<image::Image<int, 2> > itH1( *m_pGradientH1 );
	image::IteratorPosition<image::Image<int, 2> > itH2( *m_pGradientH2 );
	image::IteratorPosition<image::Image<int, 2> > itV1( *m_pGradientV1 );
	image::IteratorPosition<image::Image<int, 2> > itV2( *m_pGradientV2 );
	for (; !itH1.isAtEnd(); itH1++, itH2++) {
		int x = itH1.getIndex()[0];
		int y = itH1.getIndex()[2];
		if (((!m_stereo) || (x < width/2 - 1) || (x > width/2))
		&& ((!m_useMask) || ((*m_pMaskGradient)(x,y)))) {
			int diff = *itH1 - *itH2;
			sum += m_varH2 / (m_varH2 + diff*diff);
		}
	}
	for (; !itV1.isAtEnd(); itV1++, itV2++) {
		int x = itV1.getIndex()[0];
		int y = itV1.getIndex()[1];
		if (((!m_stereo) || (x < width/2 - 1) || (x > width/2))
		&& ((!m_useMask) || ((*m_pMaskGradient)(x,y)))) {
			int diff = *itV1 - *itV2;
			sum += m_varV2 / (m_varV2 + diff*diff);
		}
	}
	return sum / (2 * m_gradSize);
}

}
}

#endif
