/**
*  @file	OptimizerPowellBrent.h
*  @brief	Powell Optimizer with Brent line minimization
*  @author	Hari Sundar
*  @date	6/17/03
*
**/  
  
#ifndef PB_OPTIMIZER_H
#define PB_OPTIMIZER_H

#include <Optimizer.h>
#include <cmath>

static const int ITMAX = 100;
static const double GOLD = 1.618034;
static const double CGOLD = 0.3819660;
static const double ZEPS = 1.0e-10;
static const double GLIMIT = 100.0;
static const double TINY = 1.0e-20;

static void SWAP (double *a, double *b) {
  double c = *a;
  *a = *b;
  *b = c;
} 

static void SHFT (double *a, double *b, double *c, double d) {
  *a = *b;
  *b = *c;
  *c = d;
} 

static int SIGN (double x) {
  return(x != 0) ? ((x > 0) ? 1 : -1) : 0;
}


class OptimizerPowellBrent:public Optimizer 
{
  public:
    virtual bool init ();
    virtual bool optimize ();
  protected:
    double *m_dx;
    double *m_tmpx;
    double costFunction1D (double lambda) {
      for (int i = 0; i < m_numberOfParams; i++)
        m_tmpx[i] = m_parameters[i] + lambda * m_dx[i];
      m_numEvaluations++;
      return costFunction (m_tmpx);
    }
    inline void brentBracketMinimum (double *ax, double *bx, double *cx,
        double *fa, double *fb, double *fc);
    inline double brentMinimizeGivenBounds (double ax, double bx, double cx,
        double tol, double *xmin);
};

void OptimizerPowellBrent::brentBracketMinimum (double *ax, double *bx, double *cx, double *fa, double *fb, double *fc) {
  double ulim, u, r, q, fu;
  *fa = costFunction1D (*ax);
  *fb = costFunction1D (*bx);
  if (*fb > *fa) {
    SWAP (ax, bx);
    SWAP (fa, fb);
  }
  *cx = (*bx) + GOLD * (*bx - *ax);
  *fc = costFunction1D (*cx);
  while (*fb > *fc) {
    r = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    double dq = q - r;
    if (fabs (dq) < TINY)
      dq = SIGN (dq) * TINY;
    u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * dq);
    ulim = (*bx) + GLIMIT * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0) {
      fu = costFunction1D (u);
      if (fu < *fc) {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
        return;
      } else if (fu > *fb) {
        *cx = u;
        *fc = fu;
        return;
      }
      u = (*cx) + GOLD * (*cx - *bx);
      fu = costFunction1D (u);
    } else if ((*cx - u) * (u - ulim) > 0.0) {
      fu = costFunction1D (u);
      if (fu < *fc) {

        //SHFT(bx,cx,&u,*cx+GOLD*(*cx-*bx)); awf dumped -- c is useless
        SHFT (bx, cx, &u, u + GOLD * (u - *cx));
        SHFT (fb, fc, &fu, costFunction1D (u));
      }
    } else if ((u - ulim) * (ulim - *cx) >= 0.0) {
      u = ulim;
      fu = costFunction1D (u);
    } else {
      u = (*cx) + GOLD * (*cx - *bx);
      fu = costFunction1D (u);
    }
    SHFT (ax, bx, cx, u);
    SHFT (fa, fb, fc, fu);
  }
}

double OptimizerPowellBrent::brentMinimizeGivenBounds (double ax, double bx, double cx, double tol, double *xmin) {
  int iter;
  double a, b, d =
  0.0, etemp, fu, fv, fw, fx, p1, q, r, tol1, tol2, u, v, w, x, xm;
  double e = 0.0;
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = costFunction1D (x);

  //  if (verbose_) vcl_cerr << "vnl_brent f("<<x<<") \t= "<<fx <<vcl_endl;
  for (iter = 1; iter <= ITMAX; iter++) {
    xm = 0.5 * (a + b);
    tol1 = tol * fabs (x) + ZEPS;
    tol2 = 2.0 * (tol1);
    if (fabs (x - xm) <= (tol2 - 0.5 * (b - a))) {
      *xmin = x;
      return fx;
    }
    if (fabs (e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p1 = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p1 = -p1;
      q = fabs (q);
      etemp = e;
      e = d;    // Warning: The variable d has not yet been assigned a value.
      if (fabs (p1) >= fabs (0.5 * q * etemp) || p1 <= q * (a - x)
          || p1 >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));

      else {
        d = p1 / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = tol1 * SIGN (xm - x);
      }
    } else {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }
    u = (fabs (d) >= tol1 ? x + d : x + tol1 * SIGN (d));
    fu = costFunction1D (u);

    //    if (verbose_) vcl_cerr << "vnl_brent f("<<u<<") \t= "<<fu <<vcl_endl;
    if (fu <= fx) {
      if (u >= x)
        a = x;
      else
        b = x;
      SHFT (&v, &w, &x, u);
      SHFT (&fv, &fw, &fx, fu);
    } else {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  //  vcl_cerr << "Too many iterations in brent\n";
  *xmin = x;
  return fx;
}

#endif	/*  */
