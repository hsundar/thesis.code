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

class OptimizerPowellBrent:public Optimizer 
{
  public:
    virtual bool init ();
    virtual bool optimize ();
  protected:
    double *m_dx;
    double *m_tmpx;
    double costFunction1D (double lambda);
    void brentBracketMinimum (double *ax, double *bx, double *cx,
        double *fa, double *fb, double *fc);
    double brentMinimizeGivenBounds (double ax, double bx, double cx,
        double tol, double *xmin);
};

#endif	/*  */
