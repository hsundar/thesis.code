#ifndef __BSPLINE_BASIS_H_
#define __BSPLINE_BASIS_H_

class bSplineBasis {

public:
  bSplineBasis() {
    m_iOrder = 3;
    m_iNumKnots = 5;
    m_ipKnots = NULL;
    knot();
  }

  bSplineBasis(unsigned int order, unsigned int numKnots) {
    m_iOrder = order;
    m_iNumKnots = numKnots;
    m_ipKnots = NULL;
    knot();
  }

  bSplineBasis(const bSplineBasis& b) {
    m_iOrder = b.m_iOrder;
    m_iNumKnots = b.m_iNumKnots;
    m_ipKnots = NULL;
    knot();
  }
  
  bSplineBasis& operator=(const bSplineBasis &b) {
    m_iOrder = b.m_iOrder;
    m_iNumKnots = b.m_iNumKnots;
    m_ipKnots = NULL;
    knot();
		return *(this);
  }

  ~bSplineBasis() {
    if (m_ipKnots != NULL) {
      delete [] m_ipKnots;
      m_ipKnots = NULL;
    }
  };

  unsigned int getNumKnots() {
    return m_iNumKnots;
  }

  /*
  void generateBasis() {
    // generate the knots ...
    knot(m_iNumKnots, m_iOrder, m_ipKnots);
    m_dpBasis = new double*[m_ipKnots - 1];
    for (int i=0; i<m_ipKnots-1; i++) {
      m_dpBasis[i] = new double[m_ipKnots];
      basis(m_iOrder, m_iNumKnots, m_ipKnots, (double)i, m_dpBasis[i] );
    }
  }
  */

  double evaluate(double* params, double time) {
    double *bSpl = new double[m_iNumKnots];
    basis(time, bSpl);
    double sum=0.0;
    for (unsigned int i=0; i<m_iNumKnots; i++) {
      sum += params[i]*bSpl[i];
    }
    return sum;
  }

  void interpolate(double* params, int nt, double* vals) {
    double dt = 1.0/nt;
    for (int i=0; i<nt; i++) {
      vals[i] = evaluate(params, dt*i);
    }
  }

  void getBasisDot(double* lambda, int nt, double* vals) {
    double dt = 1.0/nt;
    double *params = new double[m_iNumKnots];
    for (unsigned int i=0; i<m_iNumKnots; i++) {
      params[i] = 1.0;
    }

    for (int i=0; i<nt; i++) {
      vals[i] = lambda[i]*evaluate(params, dt*i);
    }
  }

  /*
    Subroutine to generate a B-spline open knot vector with multiplicity
    equal to the order at the ends.
    */
  void knot() {  
    // m_iNumKnots, m_iOrder, m_ipKnots
    int nplusc,nplus2,i;

    nplusc = m_iNumKnots + m_iOrder;
    nplus2 = m_iNumKnots + 2;

    m_ipKnots = new int [nplusc];

    m_ipKnots[0] = 0;
    for (i = 1; i < nplusc; i++) {
      if ( (i > (-1+(int)m_iOrder)) && (i < (nplus2-1) ) )
        m_ipKnots[i] = m_ipKnots[i-1] + 1;
      else
        m_ipKnots[i] = m_ipKnots[i-1];
    }
#ifdef __DEBUG__
    std::cout << "generated knots" << std::endl;
    for (unsigned int i=0; i<m_iNumKnots+m_iOrder; i++) {
      std::cout << m_ipKnots[i] << ", ";
    }
    std::cout << std::endl;
    double *n = new double[m_iNumKnots];
    basis(0.1, n);
    std::cout << "Basis at 0.1 is ";
    for (int i=0; i<m_iNumKnots; i++) {
      std::cout << n[i] << " ";
    }
    std::cout << std::endl;
    basis(0.3, n);
    std::cout << "Basis at 0.3 is ";
    for (int i=0; i<m_iNumKnots; i++) {
      std::cout << n[i] << " ";
    }
    std::cout << std::endl;
    basis(1.4, n);
    std::cout << "Basis at 1.4 is ";
    for (int i=0; i<m_iNumKnots; i++) {
      std::cout << n[i] << " ";
    }
    std::cout << std::endl;
#endif
  }

  /*  
    Subroutine to generate B-spline basis functions for open knot vectors

    c        = order of the B-spline basis function
    d        = first term of the basis function recursion relation
    e        = second term of the basis function recursion relation
    npts     = number of defining polygon vertices
    n[]      = array containing the basis functions
               n[1] contains the basis function associated with B1 etc.
    nplusc   = constant -- npts + c -- maximum number of knot values
    t        = parameter value
    temp[]   = temporary array
    x[]      = knot vector
  */
  void basis(double tt, double *n ) {
    int c = m_iOrder;
    int npts = m_iNumKnots;
    int *x = m_ipKnots;
    int nplusc;
    int i,k;
    double d,e;

    double t = tt*(m_iNumKnots - 2);
    

    nplusc = npts + c;

    double *temp = new double[nplusc];

    // calculate the first order basis functions n[i][1]

    for (i = 0; i < nplusc-1; i++) {
      if (( t >= x[i]) && (t < x[i+1]))
        temp[i] = 1;
      else
        temp[i] = 0;
    }

    // calculate the higher order basis functions 
    for (k = 2; k <= c; k++) {
      for (i = 0; i < nplusc-k; i++) {
        if (temp[i] != 0)    /* if the lower order basis function is zero skip the calculation */
          d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
        else
          d = 0;

        if (temp[i+1] != 0)     /* if the lower order basis function is zero skip the calculation */
          e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
        else
          e = 0;

        temp[i] = d + e;
      }
    }

    if (t == (double)x[nplusc-1]) {   /*    pick up last point	*/
      temp[npts-1] = 1;
    }

    // put in n array
    // n = new float[npts];
    for (i = 0; i < npts; i++) {
      n[i] = temp[i];
    }

    // clean up
    delete [] temp;
  }

protected:
  unsigned int m_iOrder;
  unsigned int m_iNumKnots;
  int* m_ipKnots;
};

#endif
