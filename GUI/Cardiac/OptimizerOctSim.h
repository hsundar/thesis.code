#ifndef __OPTIMIZER_OCT_SIM_H__
#define __OPTIMIZER_OCT_SIM_H__

#include "oct_sim.h"
#include "TransMatrix.h"
#include "OptimizerGradientDescent.h"
#include "OptimizerPowellBrent.h"

class OptimizerOctGD : public OptimizerGradientDescent
{
  public:
    virtual double costFunction (double *tmpParams)
    {
      // std::cout << "Evaluating cost function" << std::endl;
      // 6 params are -> tx, ty, tz, rx, ry, rz

      // convert params into transform and give it to oct_sim
      TransMatrix trans = TransMatrix::TransRot(tmpParams);
      
      m_costFunctionValue = -(sm->getMI(trans.GetDataPointer()));
      // verbose();
      return m_costFunctionValue; 
    }

    void setSimilarityMeasure(oct_sim* sim) { 
      sm = sim; 
    }

    virtual void verbose ()
    {
      printf ("[ ");
      for (int i = 0; i < m_numberOfParams; i++)
        printf ("%f ", m_parameters[i]);
      printf ("] -> %f\n", m_costFunctionValue);
    } 
  
  private:
    oct_sim *sm;
};

class OptimizerOctPB : public OptimizerPowellBrent
{
  public:
    virtual double costFunction (double *tmpParams)
    {  
      //std::cout << "In cost function: "; // << std::endl;
      // 6 params are -> tx, ty, tz, rx, ry, rz
      // convert params into transform and give it to oct_sim
      TransMatrix trans = TransMatrix::TransRot(tmpParams);
      
      // sm.CalculateMeasures ();

      double cost = sm->getMI(trans.GetDataPointer());
      // std::cout << cost << std::endl;
      return -cost;
    }
    virtual void verbose ()
    {
      printf ("[ ");
      for (int i = 0; i < m_numberOfParams; i++)
        printf ("%f ", m_parameters[i]);
      printf ("] -> %f\n", m_costFunctionValue);
    } 
    void setSimilarityMeasure(oct_sim* sim) { 
      sm = sim; 
    }
  private:
    oct_sim *sm;
};


#endif

