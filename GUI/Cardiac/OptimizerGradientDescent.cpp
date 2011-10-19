/**
 *  @file      OptimizerGradientDescent.cpp
 *  @author    Hari Sundar
 *  @date      06/06/03 
 * 
 **/  
#include "OptimizerGradientDescent.h"
 
#include <iostream>

bool OptimizerGradientDescent::init ()
{
  // (re)initialize memory for work parameters and m_workResults
  if (m_workParams)
    delete[]m_workParams;

  m_workParams = new double[m_numberOfParams];

  if (m_workResults)
    delete[]m_workResults;

  m_workResults = new double[m_numberOfParams];

  // set default options
  m_maxIterations = 10000;

  m_maxEvaluations = 1000000;

  m_functionTolerance = 0.000001;

  m_paramsTolerance = 0.000001;

  m_learningRate = 1.0;

  return true;

}


bool OptimizerGradientDescent::optimize ()
{
  std::cout << "In Optimize" << std::endl;

  m_isOptimizing = true;

  int dims = m_numberOfParams;

  m_numIterations = 0;

  m_numEvaluations = 0;

  bool downscaled = false, better = false;

  double lastvalue = 0, currvalue = 0;

  m_scale = 1.0;


  std::cout << "In Optimize: Initialize " << m_numberOfParams << std::endl;
  // initialize work parameter space
  for (int i = 0; i < m_numberOfParams; i++)
    m_workParams[i] = m_parameters[i];


  std::cout << "In Optimize: Calling costFx" << std::endl;
  // evaluate the current position first
  double c = costFunction (m_workParams);

  m_costFunctionValue = c;

  // negate result if we are minimizing
  currvalue = (m_minimize) ? -c : c;


  // main optimization loop
  while ((m_numIterations < m_maxIterations)
      && (m_numEvaluations <
        m_maxEvaluations) 
      &&((m_numIterations == 0)
        || (downscaled)
        || (currvalue - lastvalue >
          m_functionTolerance)))
  {


    // save last result
    lastvalue = currvalue;

    //int bestpos = 0;


    // calculate gradient
    for (int i = 0; i < dims; i++)
    {


      // move into specific direction
      m_workParams[i] += m_stepsizes[i] * m_scale;


      // calculate cost function
      double c = costFunction (m_workParams);

      // negate result if we are minimizing
      if (m_minimize)
        c = -c;

      // set gradient entry
      m_workResults[i] =
        (c - lastvalue) / (m_stepsizes[i] * m_scale);


      // reset work parameter
      m_workParams[i] = m_parameters[i];


      m_numEvaluations++;

    }

    // apply the step
    for (int i = 0; i < dims; i++)
      m_workParams[i] += m_workResults[i] * m_scale * m_learningRate;

    // evaluate the new position
    double c = costFunction (m_workParams);

    // negate result if we are minimizing
    if (m_minimize)
      c = -c;

    if (c > currvalue)
    {

      for (int i = 0; i < dims; i++)
        m_parameters[i] = m_workParams[i];

      m_costFunctionValue = (m_minimize) ? -c : c;

      currvalue = c;

      better = true;

      downscaled = false;

    }
    else
    {

      m_scale *= 0.5;

      downscaled = true;

      // for now, use the scaling factor as parameter abort criterion
      if (m_scale < m_paramsTolerance)
        break;

    }


    m_numIterations++;

    verbose ();

  }


  m_isOptimizing = false;

  return better;

}


