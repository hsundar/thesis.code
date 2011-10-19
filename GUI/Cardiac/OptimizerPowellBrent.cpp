/**
 *  @file	OptimizerPowellBrent.cpp
 *  @brief	Implements the OptimizerPowellBrent class
 *  @author	Hari Sundar
 *  @date	6/17/03
 *
 **/  

#include "OptimizerPowellBrent.h"

bool OptimizerPowellBrent::init () {

  // (re)initialize memory for work parameters and results
  if (m_workParams)
    delete[]m_workParams;
  m_workParams = new double[m_numberOfParams];
  if (m_workResults)
    delete[]m_workResults;
  m_workResults = new double[2 * m_numberOfParams];

  // set default options
  m_maxIterations = 100;
  m_maxEvaluations = 10000;
  m_functionTolerance = 0.001;
  m_paramsTolerance = 0.001;
  return true;
}

bool OptimizerPowellBrent::optimize () {
  m_isOptimizing = true;
  int i, j;
  int n = m_numberOfParams;
  m_numEvaluations = 0;
  m_numIterations = 0;

  // initialize parameter space
  for (i = 0; i < m_numberOfParams; i++)
    m_workParams[i] = this->m_parameters[i];
  double *p = m_parameters;
  double *pt = m_workParams;
  double *ptt = new double[n];
  double *xit = new double[n];
  m_dx = xit;
  m_tmpx = new double[n];
  double linmin_ftol = 1e-4;

  // initialize n * n array and set the identity matrix

  // hs::array2d < double >xi (n, n);
  double *xi = new double[n*n];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      xi[n*i+j] = i == j ? 1.0 : 0.0;

  double fret = costFunction (p);
  m_costFunctionValue = fret;
  m_numEvaluations++;

  for (;;) {
    double fp = fret;
    int ibig = 0;
    double del = 0.0;
    for (i = 0; i < n; i++) {
      // xit = ith column of xi
      for (j = 0; j < n; ++j)
        xit[j] = xi[n*j+i];
      double fptt = fret;

      // 1D minimization along xi
      double ax = 0.0;
      double xx = 1.0;
      double bx, fa, fb, fc;
      brentBracketMinimum (&ax, &xx, &bx, &fa, &fb, &fc);
      fret = brentMinimizeGivenBounds (bx, xx, ax, linmin_ftol, &xx);

      // set the new position
      for (j = 0; j < n; j++)
        p[j] = m_parameters[j] + xx * m_dx[j];

      // Now p is minimizer along xi
      if (fabs (fptt - fret) > del) {
        del = fabs (fptt - fret);
        ibig = i;
      }
    }
    if (2.0 * fabs (fp - fret) <=
        m_functionTolerance * (fabs (fp) + fabs (fret))) {

      //                      vnl_matlab_print(vcl_cerr, xi, "xi");
      //                      return CONVERGED_FTOL;
      break;
    }
    if (m_numEvaluations >= m_maxEvaluations) {

      //                      return FAILED_TOO_MANY_ITERATIONS;
      break;
    }
    for (int j = 0; j < n; ++j) {
      ptt[j] = 2.0 * p[j] - pt[j];
      xit[j] = p[j] - pt[j];
      pt[j] = p[j];
    } double fptt = costFunction (ptt);
    m_costFunctionValue = fret;
    m_numEvaluations++;
    if (fptt < fp) {
      double t =
      2.0 * (fp - 2.0 * fret +
             fptt) * ((fp - fret - del) * (fp - fret - del)) -
      del * ((fp - fptt) * (fp - fptt));
      if (t < 0.0) {
        double ax = 0.0;
        double xx = 1.0;
        double bx, fa, fb, fc;
        brentBracketMinimum (&ax, &xx, &bx, &fa, &fb, &fc);
        fret =
        brentMinimizeGivenBounds (bx, xx, ax, linmin_ftol, &xx);

        // set the new position
        for (int j = 0; j < n; j++)
          p[j] = m_parameters[j] + xx * m_dx[j];
        for (int j = 0; j < n; j++) {
          xi[n*j + ibig] = xi[n*j + n - 1];
          xi[n*j + n - 1] = xit[j];
        }
      }
    }
    m_numIterations++;
    verbose ();
  } // done infinite loop  
  delete [] ptt;
  delete [] xit;
  delete [] xi;
  delete [] m_tmpx;
  m_isOptimizing = false;
  return true;
}


