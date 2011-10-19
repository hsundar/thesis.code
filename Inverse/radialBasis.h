#ifndef __RADIAL_BASIS_H_
#define __RADIAL_BASIS_H_

#include "Point.h"

class radialBasis {
public:
  radialBasis() {
    m_ptCenter = Point();
    m_ptSigmaSq = Point(1.0, 1.0, 1.0);
  }

  radialBasis(Point ctr, Point sig2) {
    m_ptCenter = ctr;
    m_ptSigmaSq = sig2;
  }

  ~radialBasis() {};

	radialBasis(const radialBasis& other) {
		m_ptCenter = other.m_ptCenter;
		m_ptSigmaSq = other.m_ptSigmaSq;
	}

  double getValue(Point x) {
    // return the value at point x ...
    Point d = x - m_ptCenter;
    double val = d.x()*d.x()/m_ptSigmaSq.x() + d.y()*d.y()/m_ptSigmaSq.y() + d.z()*d.z()/m_ptSigmaSq.z();
    val *= -0.5;
    return exp(val);
  }

  void getCenter(double &x, double &y, double& z) {
    x = m_ptCenter.x(); y = m_ptCenter.y(); z = m_ptCenter.z();
  }

protected:
  Point m_ptCenter;
  Point m_ptSigmaSq;

};

#endif
