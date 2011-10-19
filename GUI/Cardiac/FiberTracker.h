#ifndef __FIBER_TRACKER_H_
#define __FIBER_TRACKER_H_

#include <vector>

#include "Field.h"
#include "Point.h"

#define FA_THRESHOLD 0.1 
#define max_iter 500000
#define step_sz 1.0

class SoSeparator;

class FiberTracker {
  public:
    enum Color { Gray, Normal, Chiral};
    FiberTracker ();
    virtual ~FiberTracker();

    void setFractionalAnisotropy(Field<float, 1> *fa) { m_pFA = fa; };
    void setPrincipalDirection(Field<float, 3> *pd) { m_pPD = pd; };

    void setSeedDensity (float rho) { m_fSeedDensity = rho; };
    void setColorFunction(Color clr=Gray) {m_iColor = clr;};

    SoSeparator* track();

  protected:
    Field<float, 1>	*m_pFA;
    Field<float, 3>	*m_pPD;

    float		m_fSeedDensity;

    int			m_iColor;
    // Containers ...
    std::vector<Point> seeds;
    std::vector<Point> fiber;
    std::vector<Point> norms;

  private:
    Point getVelocity(int i, int j, int k);
    Point getVelocity(Point p, Point lastV);
    Point getVelocityInterpolated(Point p);
    unsigned int pt2Col(Point p, Point q);
    unsigned int float2Col(float p);

};

#endif
