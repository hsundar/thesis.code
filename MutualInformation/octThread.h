#ifndef __OCT_THREAD_H_
#define __OCT_THREAD_H_

#include <QThread>
#include <QMutex>
#include <QWaitCondition>

#include <vector>

#include "Volume.h"

/*****************************************************************************
 * Computes the similarity between two images, subject to an affine transform
 * between them. The Class is threaded, and therefore can work in parallel on
 * multi-core machines.
 *
 * The parent thread will define domains on which the individual threads will
 * operate.
 *
 *****************************************************************************/


class octThread : public QThread {
     Q_OBJECT

 public:
     
     enum similarityMeasureType { SSD, NCC, MI, NMI };  
       
     octThread(QObject *parent = 0);
     ~octThread();
     
     void setSourceImage(Volume * img) { m_sourceImg = img; }
     void setTargetImage(Volume * img) { m_targetImg = img; }

     void setDomain(unsigned int x1, unsigned int y1, unsigned int z1, unsigned int x2, unsigned int y2, unsigned int z2) {
       m_iMinX = x1; m_iMinY = y1; m_iMinZ = z1;
       m_iMaxX = x2; m_iMaxY = y2; m_iMaxZ = z2;
     }

     void setDomain (unsigned int o1, unsigned int o2) {
        m_iMinOct = o1;
        m_iMaxOct = o2;
     }

     void setSimilarityMeasure(int m) { m_iSimilarityType = m; }
     void setHistBits(unsigned int h) { m_histBits = h;}
     void useOctree(bool flag) { m_bUseOctree = flag; }

     void compute(double *trans);

  protected:
     void initHists(); // initilizes the histograms. Required for MI, NMI
     
     double getSSD(double *transMat);
     double getNCC(double *transMat);
     double getNMI(double *transMat);
     double getMI (double *transMat);

 protected:
     void run();

 private:
     // Thread related ...
     QMutex mutex;
     QWaitCondition condition;
     
     bool abort;

     // oct sim related ...
     Volume *		m_sourceImg;
     Volume *		m_targetImg;

     // The histograms ...
     std::vector<unsigned int>	m_sourceHist;
     std::vector<unsigned int>	m_targetHist;

     unsigned int*		m_jointHist;
     unsigned int  		m_histBits;

     double*			m_transform;

     int 			m_iSimilarityType;
     bool			m_bUseOctree;

     unsigned int		m_iMinX;
     unsigned int		m_iMinY;
     unsigned int		m_iMinZ;
     unsigned int		m_iMinOct;
     
     unsigned int		m_iMaxX;
     unsigned int		m_iMaxY;
     unsigned int		m_iMaxZ;
     unsigned int		m_iMaxOct;
};

#endif

