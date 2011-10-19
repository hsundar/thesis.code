#include "oct_sim.h"
#include <iostream>

oct_sim::oct_sim() {
  m_sourceImg = NULL;
  m_targetImg = NULL;

  m_histBits = 6; // 64 bins ...
  m_jointHist = NULL;

  m_bUseOctree = false;
}

oct_sim::~oct_sim() {
  if (m_jointHist != NULL) {
    delete [] m_jointHist;
    m_jointHist = NULL;
  }
  m_sourceHist.clear();
  m_targetHist.clear();
}

void oct_sim::initHists() {
  m_sourceHist.clear();
  m_targetHist.clear();

  unsigned int histSize = 1 << m_histBits;

  m_sourceHist.resize(histSize);
  m_targetHist.resize(histSize);

  m_jointHist = new unsigned int [histSize*histSize];
}

double oct_sim::getSSD(double *transMat) {
  double sum = 0.0;
  if (!m_bUseOctree) {
    for (int k=0; k< m_targetImg->GetSize(2); k++) {
      for (int j=0; j< m_targetImg->GetSize(1); j++) {
        for (int i=0; i< m_targetImg->GetSize(0); i++) {
          Point _tar(i,j,k);
          Point _src = Point::TransMatMultiply(transMat, _tar);
          double x = m_sourceImg->GetAt(_src) - m_targetImg->GetAt(_tar);
          sum +=  x*x;
        } 
      } 
    } 
  } else {
    // loop over all elements ...
    for (unsigned int i=0; i< m_targetImg->octSize(); i++) {
      Point _tar = m_targetImg->getNodeCenter(i);
      Point _src = Point::TransMatMultiply(transMat, _tar);
      double x = m_sourceImg->GetAt(_src) - m_targetImg->GetAt(_tar);
      sum +=  x*x;
    }
  }

  return sum;
}

double oct_sim::getNCC(double *transMat) {
  double sum1=0., sum2=0., sum3=0., sum4=0., sum5=0.;
  int size = 0;
  if (!m_bUseOctree) {
    for (int k=0; k< m_targetImg->GetSize(2); k++) {
      for (int j=0; j< m_targetImg->GetSize(1); j++) {
        for (int i=0; i< m_targetImg->GetSize(0); i++) {
          Point _tar(i,j,k);
          Point _src = Point::TransMatMultiply(transMat, _tar);
          double pixel1 = m_sourceImg->GetAt(_src);
          double pixel2 = m_targetImg->GetAt(_tar);

          sum1 += pixel1;
          sum2 += pixel2;
          sum3 += pixel1 * pixel1;
          sum4 += pixel2 * pixel2;
          sum5 += pixel1 * pixel2;
          size++;
        } 
      } 
    }
  } else {
    for (unsigned int i=0; i< m_targetImg->octSize(); i++) {
      Point _tar = m_targetImg->getNodeCenter(i);
      Point _src = Point::TransMatMultiply(transMat, _tar);

      double pixel1 = m_sourceImg->GetAt(_src);
      double pixel2 = m_targetImg->GetAt(_tar);

      sum1 += pixel1;
      sum2 += pixel2;
      sum3 += pixel1 * pixel1;
      sum4 += pixel2 * pixel2;
      sum5 += pixel1 * pixel2;
      size++;
    }
  }


  double var1 = 0, var2 = 0, corr = 0;

  corr = ((double)sum5) - ((double)(sum1 * sum2)) / (double)size;
  var1 = ((double)sum3) - ((double)(sum1 * sum1)) / (double)size;
  var2 = ((double)sum4) - ((double)(sum2 * sum2)) / (double)size;

  return corr / (sqrt(var1) * sqrt(var2));
}

double oct_sim::getNMI(double *transMat) {
  int histSize = 1 << m_histBits;

  // clear the histograms ...
  for (int i=0; i<histSize; i++) {
    m_sourceHist[i] = 0;
    m_targetHist[i] = 0;
    for (int j=0; j<histSize; j++) 
      m_jointHist[i*histSize + j] = 0;
  }

  // first compute the histograms ...
  unsigned int bitsUsed = m_sourceImg->GetVolBitsStored();
  int size = 0;

  if (!m_bUseOctree) {
    for (int k=0; k< m_targetImg->GetSize(2); k++) {
      for (int j=0; j< m_targetImg->GetSize(1); j++) {
        for (int i=0; i< m_targetImg->GetSize(0); i++) {
          Point _tar(i,j,k);
          Point _src = Point::TransMatMultiply(transMat, _tar);
          unsigned int pixel1 = (unsigned int)m_sourceImg->GetAt(_src);
          unsigned int pixel2 = (unsigned int)m_targetImg->GetAt(_tar);

          // scale values for histogram indexing
          unsigned int histIndex1;
          unsigned int histIndex2;

          histIndex1 = pixel1 >> (bitsUsed - m_histBits);
          histIndex2 = pixel2 >> (bitsUsed - m_histBits);

          m_jointHist[histIndex2*histSize + histIndex2]++;

          m_sourceHist[histIndex1]++;
          m_targetHist[histIndex2]++;

          size++;
        } 
      } 
    }
  } else {
    for (unsigned int i=0; i< m_targetImg->octSize(); i++) {
      Point _tar = m_targetImg->getNodeCenter(i);
      Point _src = Point::TransMatMultiply(transMat, _tar);

      unsigned int pixel1 = (unsigned int)m_sourceImg->GetAt(_src);
      unsigned int pixel2 = (unsigned int)m_targetImg->GetAt(_tar);

      // scale values for histogram indexing
      unsigned int histIndex1;
      unsigned int histIndex2;

      histIndex1 = pixel1 >> (bitsUsed - m_histBits);
      histIndex2 = pixel2 >> (bitsUsed - m_histBits);

      m_jointHist[histIndex2*histSize + histIndex2]++;

      m_sourceHist[histIndex1]++;
      m_targetHist[histIndex2]++;

      size++;
    }
  }

  double p;
  double src_entropy=0., tar_entropy=0., joint_entropy=0.;
  for (int y = 0; y < histSize; y++) {
    p = ((double) m_sourceHist[y]) / ((double) size);
    if (p > 0.0) src_entropy -= p * log(p);
    p = ((double) m_targetHist[y]) / ((double) size);
    if (p > 0.0) tar_entropy -= p * log(p);
    for (int x = 0; x < histSize; x++) {
      p = ((double) m_jointHist[y*histSize + x]) / ((double) size);
      if (p > 0.0) joint_entropy -= p * log(p);
    }
  } // for (int y = 0; y < histSize; y++)



  // std::cout << "Entropies " << joint_entropy << " " << src_entropy << " " << tar_entropy << std::endl;
  if ( (src_entropy + tar_entropy) > 0.0)
    return  (2. - 2. * joint_entropy / (src_entropy + tar_entropy) );
  else
    return 0.0;
}


double oct_sim::getMI (double *transMat) {
  int histSize = 1 << m_histBits;

  // clear the histograms ...
  for (int i=0; i<histSize; i++) {
    m_sourceHist[i] = 0;
    m_targetHist[i] = 0;
    for (int j=0; j<histSize; j++) 
      m_jointHist[i*histSize + j] = 0;
  }
  // first compute the histograms ...
  unsigned int bitsUsed = m_sourceImg->GetVolBitsStored();
  int size = 0;

  if (!m_bUseOctree) {
    for (int k=0; k< m_targetImg->GetSize(2); k++) {
      for (int j=0; j< m_targetImg->GetSize(1); j++) {
        for (int i=0; i< m_targetImg->GetSize(0); i++) {
          Point _tar(i,j,k);
          Point _src = Point::TransMatMultiply(transMat, _tar);
          unsigned int pixel1 = (unsigned int)m_sourceImg->GetAt(_src);
          unsigned int pixel2 = (unsigned int)m_targetImg->GetAt(_tar);

          // scale values for histogram indexing
          unsigned int histIndex1;
          unsigned int histIndex2;

          histIndex1 = pixel1 >> (bitsUsed - m_histBits);
          histIndex2 = pixel2 >> (bitsUsed - m_histBits);

          m_jointHist[histIndex2*histSize + histIndex2]++;

          m_sourceHist[histIndex1]++;
          m_targetHist[histIndex2]++;

          size++;
        } 
      } 
    }
  } else {
    for (unsigned int i=0; i< m_targetImg->octSize(); i++) {
      Point _tar = m_targetImg->getNodeCenter(i);
      //std::cout << "Target is " << _tar.x() << ", " << _tar.y() << ", " << _tar.z() << std::endl;
      Point _src = Point::TransMatMultiply(transMat, _tar);
      //std::cout << "Source is " << _src.x() << ", " << _src.y() << ", " << _src.z() << std::endl << std::endl;

      unsigned int pixel1 = (unsigned int)m_sourceImg->GetAt(_src);
      unsigned int pixel2 = (unsigned int)m_targetImg->GetAt(_tar);

      // scale values for histogram indexing
      unsigned int histIndex1;
      unsigned int histIndex2;

      histIndex1 = pixel1 >> (bitsUsed - m_histBits);
      histIndex2 = pixel2 >> (bitsUsed - m_histBits);

      m_jointHist[histIndex2*histSize + histIndex2]++;

      m_sourceHist[histIndex1]++;
      m_targetHist[histIndex2]++;

      size++;
    }
  }

  double p;
  double src_entropy=0., tar_entropy=0., joint_entropy=0.;
  for (int y = 0; y < histSize; y++) {
    p = ((double) m_sourceHist[y]) / ((double) size);
    if (p > 0.0) src_entropy -= p * log(p);
    p = ((double) m_targetHist[y]) / ((double) size);
    if (p > 0.0) tar_entropy -= p * log(p);
    for (int x = 0; x < histSize; x++) {
      p = ((double) m_jointHist[y*histSize + x]) / ((double) size);
      if (p > 0.0) joint_entropy -= p * log(p);
    }
  }
  /*
     if (m_normalizeMI) {
     if (src_entropy + tar_entropy > 0.0)
     return  (2 - 2 * joint_entropy / (src_entropy + tar_entropy) );
     }
     else
     */
  // std::cout << " " << src_entropy << " " << tar_entropy << " " << joint_entropy << std::endl; 
  return  src_entropy + tar_entropy - joint_entropy;
}
