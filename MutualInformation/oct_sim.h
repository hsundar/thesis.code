#ifndef __OCT_SIM_H_
#define __OCT_SIM_H_

// simple similarity measure class.

// currently computes, SSD, NCC, MI and NMI

#include <vector>
#include "Volume.h"

class oct_sim {
  Volume *		m_sourceImg;
  Volume *		m_targetImg;

  // The histograms ...
  std::vector<unsigned int>	m_sourceHist;
  std::vector<unsigned int>	m_targetHist;

  unsigned int*			m_jointHist;
  unsigned int  		m_histBits;
  
  bool	m_bUseOctree;
  public:
    oct_sim();
    virtual ~oct_sim();

    void setSourceImage(Volume * img) { m_sourceImg = img; }
    void setTargetImage(Volume * img) { m_targetImg = img; }

    void setHistBits(unsigned int h) { m_histBits = h;}

    void initHists(); // initilizes the histograms. Required for MI, NMI
    void useOctree(bool flag) { m_bUseOctree = flag; }

    double getSSD(double *transMat);
    double getNCC(double *transMat);
    double getNMI(double *transMat);
    double getMI (double *transMat);
};

#endif
