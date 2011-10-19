/**
*  @file	OptimizerGradientDescent.h
*  @brief	Gradient Descent Optimizer
*  @author	Hari Sundar
*  @date	06/06/03
*
**/  
  
#ifndef OPTIMIZERGRADIENTDESCENT_H
#define OPTIMIZERGRADIENTDESCENT_H
  
#include <Optimizer.h>
  
class OptimizerGradientDescent:public Optimizer 
{
  public:
    virtual bool init ();
    virtual bool optimize ();
    double getLearningRate ()
    {
      return m_learningRate;
    };
    void setLearningRate (double val)
    {
      m_learningRate = val;
    };
  protected:
    double m_learningRate;
};

#endif	/*  */
