/**
 *  @file	OptimizerPowellBrent.cpp
 *  @brief	Implements the OptimizerPowellBrent class
 *  @author	Hari Sundar
 *  @date	6/17/03
 *
 **/  

#include "OptimizerPowellBrent.h"
#include <cmath>

static const int ITMAX = 100;
static const double GOLD = 1.618034;
static const double CGOLD = 0.3819660;
static const double ZEPS = 1.0e-10;
static const double GLIMIT = 100.0;
static const double TINY = 1.0e-20;
static void SWAP (double *a, double *b)
{
  double c = *a;
  *a = *b;
  *b = c;
} static void SHFT (double *a, double *b, double *c, double d)
{
  *a = *b;
  *b = *c;
  *c = d;
} static int SIGN (double x)
{
  return (x != 0) ? ((x > 0) ? 1 : -1) : 0;
}
bool OptimizerPowellBrent::init ()
{

  // (re)initialize memory for work parameters and results
  if (m_workParams)
    delete[]m_workParams;
  m_workParams = new double[m_numberOfParams];
  if (m_workResults)
    delete[]m_workResults;
  m_workResults = new double[2 * m_numberOfParams];

  // set default options
  m_maxIterations = 10000;
  m_maxEvaluations = 1000000;
  m_functionTolerance = 0.000001;
  m_paramsTolerance = 0.000001;
  return true;
}
double OptimizerPowellBrent::costFunction1D (double lambda)
{
  for (int i = 0; i < m_numberOfParams; i++)
    m_tmpx[i] = m_parameters[i] + lambda * m_dx[i];
  m_numEvaluations++;
  return costFunction (m_tmpx);
}
  void OptimizerPowellBrent::brentBracketMinimum 
(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc)
{
  double ulim, u, r, q, fu;
  *fa = costFunction1D (*ax);
  *fb = costFunction1D (*bx);
  if (*fb > *fa)
  {
    SWAP (ax, bx);
    SWAP (fa, fb);
  }
  *cx = (*bx) + GOLD * (*bx - *ax);
  *fc = costFunction1D (*cx);
  while (*fb > *fc)
  {
    r = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    double dq = q - r;
    if (fabs (dq) < TINY)
      dq = SIGN (dq) * TINY;
    u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * dq);
    ulim = (*bx) + GLIMIT * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0)
    {
      fu = costFunction1D (u);
      if (fu < *fc)
      {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
        return;
      }
      else if (fu > *fb)
      {
        *cx = u;
        *fc = fu;
        return;
      }
      u = (*cx) + GOLD * (*cx - *bx);
      fu = costFunction1D (u);
    }
    else if ((*cx - u) * (u - ulim) > 0.0)
    {
      fu = costFunction1D (u);
      if (fu < *fc)
      {

        //SHFT(bx,cx,&u,*cx+GOLD*(*cx-*bx)); awf dumped -- c is useless
        SHFT (bx, cx, &u, u + GOLD * (u - *cx));
        SHFT (fb, fc, &fu, costFunction1D (u));
      }
    }
    else if ((u - ulim) * (ulim - *cx) >= 0.0)
    {
      u = ulim;
      fu = costFunction1D (u);
    }
    else
    {
      u = (*cx) + GOLD * (*cx - *bx);
      fu = costFunction1D (u);
    }
    SHFT (ax, bx, cx, u);
    SHFT (fa, fb, fc, fu);
  }
}
  double OptimizerPowellBrent::brentMinimizeGivenBounds 
(double ax, double bx, double cx, double tol, double *xmin)
{
  int iter;
  double a, b, d =
    0.0, etemp, fu, fv, fw, fx, p1, q, r, tol1, tol2, u, v, w, x, xm;
  double e = 0.0;
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = costFunction1D (x);

  //  if (verbose_) vcl_cerr << "vnl_brent f("<<x<<") \t= "<<fx <<vcl_endl;
  for (iter = 1; iter <= ITMAX; iter++)
  {
    xm = 0.5 * (a + b);
    tol1 = tol * fabs (x) + ZEPS;
    tol2 = 2.0 * (tol1);
    if (fabs (x - xm) <= (tol2 - 0.5 * (b - a)))
    {
      *xmin = x;
      return fx;
    }
    if (fabs (e) > tol1)
    {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p1 = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p1 = -p1;
      q = fabs (q);
      etemp = e;
      e = d;		// Warning: The variable d has not yet been assigned a value.
      if (fabs (p1) >= fabs (0.5 * q * etemp) || p1 <= q * (a - x)
          || p1 >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));

      else
      {
        d = p1 / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = tol1 * SIGN (xm - x);
      }
    }
    else
    {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }
    u = (fabs (d) >= tol1 ? x + d : x + tol1 * SIGN (d));
    fu = costFunction1D (u);

    //    if (verbose_) vcl_cerr << "vnl_brent f("<<u<<") \t= "<<fu <<vcl_endl;
    if (fu <= fx)
    {
      if (u >= x)
        a = x;
      else
        b = x;
      SHFT (&v, &w, &x, u);
      SHFT (&fv, &fw, &fx, fu);
    }
    else
    {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x)
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
        v = u;
        fv = fu;
      }
    }
  }

  //  vcl_cerr << "Too many iterations in brent\n";
  *xmin = x;
  return fx;
}
bool OptimizerPowellBrent::optimize ()
{
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

  for (;;)
  {
    double fp = fret;
    int ibig = 0;
    double del = 0.0;
    for (i = 0; i < n; i++)
    {
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
      if (fabs (fptt - fret) > del)
      {
        del = fabs (fptt - fret);
        ibig = i;
      }
    }
    if (2.0 * fabs (fp - fret) <=
        m_functionTolerance * (fabs (fp) + fabs (fret)))
    {

      //                      vnl_matlab_print(vcl_cerr, xi, "xi");
      //                      return CONVERGED_FTOL;
      break;
    }
    if (m_numEvaluations >= m_maxEvaluations)
    {

      //                      return FAILED_TOO_MANY_ITERATIONS;
      break;
    }
    for (int j = 0; j < n; ++j)
    {
      ptt[j] = 2.0 * p[j] - pt[j];
      xit[j] = p[j] - pt[j];
      pt[j] = p[j];
    } double fptt = costFunction (ptt);
    m_costFunctionValue = fret;
    m_numEvaluations++;
    if (fptt < fp)
    {
      double t =
        2.0 * (fp - 2.0 * fret +
            fptt) * ((fp - fret - del) * (fp - fret - del)) -
        del * ((fp - fptt) * (fp - fptt));
      if (t < 0.0)
      {
        double ax = 0.0;
        double xx = 1.0;
        double bx, fa, fb, fc;
        brentBracketMinimum (&ax, &xx, &bx, &fa, &fb, &fc);
        fret =
          brentMinimizeGivenBounds (bx, xx, ax, linmin_ftol, &xx);

        // set the new position
        for (int j = 0; j < n; j++)
          p[j] = m_parameters[j] + xx * m_dx[j];
        for (int j = 0; j < n; j++)
        {
          xi[n*j + ibig] = xi[n*j + n - 1];
          xi[n*j + n - 1] = xit[j];
        } }
    }
    m_numIterations++;
    verbose ();
  } 
  delete [] ptt;
  delete [] xit;
  delete [] xi;
  delete [] m_tmpx;
  m_isOptimizing = false;
  return true;
}


