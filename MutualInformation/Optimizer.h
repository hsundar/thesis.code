/**
*  @file	Optimizer.h
*  @brief	The base class for Optimizers
*  @author	Hari Sundar
*  @date	5/29/03
*
*  This is supposed to be a base class for optimizers as such. However, because of limited use currently
*  it is being written as a non-linear optimizer. It can be later on rewritten as a non-linear optimizer 
*  derived from a more generic Optimizer class.
*
**/  
  
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

class Optimizer 
{
  public:
    /** @name Constructors and Destructors **/ 
    //@{
    Optimizer ();
    virtual ~ Optimizer ();

    //@}
    virtual bool init () = 0;
    virtual double costFunction (double *tmpParams) = 0;
    virtual bool optimize () = 0;
    virtual void verbose () = 0;

    /** @name access methods **/ 
    //@{
    void setParamsTolerance (double x)
    {
      m_paramsTolerance = x;
    } void setFunctionTolerance (double x)
    {
      m_functionTolerance = x;
    } void setNumberOfIterations (int n)
    {
      m_maxIterations = n;
    } void setNumberOfEvaluations (int n)
    {
      m_maxEvaluations = n;
    } void setMinimize (bool flag)
    {
      m_minimize = flag;
    } double getParamsTolerance () const
    {
      return m_paramsTolerance;
    }
    double getFunctionTolerance () const
    {
      return m_functionTolerance;
    }
    int getMaximumNumberOfIterations () const
    {
      return m_maxIterations;
    }
    int getMaximumNumberOfEvaluations () const
    {
      return m_maxEvaluations;
    }
    int getFinalNumberOfIterations () const
    {
      return m_numIterations;
    }
    int getFinalNumberOfEvaluations () const
    {
      return m_numEvaluations;
    }
    double getCostFunctionValue () const
    {
      return m_costFunctionValue;
    };
    void setInitialParameters (double *parameters, double *stepsizes,
        int dimensions);
    double *getParameters () const
    {
      return m_parameters;
    };
    bool isOptimizing ()const
    {
      return m_isOptimizing;
    };

    //@}
  protected:
    /** @name Data members **/ 
    //@{
    double m_costFunctionValue;
    double m_paramsTolerance;
    double m_functionTolerance;
    int m_maxIterations;
    int m_numIterations;
    int m_maxEvaluations;
    int m_numEvaluations;
    double *m_parameters;
    double *m_stepsizes;
    double *m_workParams;
    double *m_workResults;
    double m_scale;
    int m_numberOfParams;
    bool m_isOptimizing;
    bool m_minimize;

    //@}
};
#endif	/*  */

