/**
*  @file	BestNeighbor.cpp
*  @brief	Implements the OptimizerBestNeighbor class
*  @author	Hari Sundar
*  @date	6/2/03
*
**/  
  
#include <OptimizerBestNeighbor.h>

bool OptimizerBestNeighbor::init ()
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
bool OptimizerBestNeighbor::optimize ()
{
  m_isOptimizing = true;
  int dims = m_numberOfParams;
  m_numEvaluations = 0;
  m_numIterations = 0;
  bool downscaled = false, better = false;
  double lastvalue = 0, currvalue = 0;
  m_scale = 1;

  // initialize parameter space
  for (int i = 0; i < m_numberOfParams; i++)
    m_workParams[i] = this->m_parameters[i];

  // evaluate the current position first
  double c = costFunction (m_workParams);
  m_costFunctionValue = c;

  // negate result if we are minimizing
  currvalue = (m_minimize) ? -c : c;

  // main optimization loop
  while ((m_numIterations < m_maxIterations)
      && (m_numEvaluations <
        m_maxEvaluations)  &&((m_numIterations == 0)
          || (downscaled)
          || (currvalue - lastvalue >
            m_functionTolerance)))
  {

    // save last result
    lastvalue = currvalue;
    int bestpos = 0;

    // evaluate neighborhood
    for (int i = 0; i < 2 * dims; i++)
    {

      // move into specific direction
      m_workParams[i % dims] +=
        m_stepsizes[i % dims] * m_scale * (i < dims ? 1.0 : -1.0);

      // calculate cost function
      double c = costFunction (m_workParams);

      // negate result if we are minimizing
      m_workResults[i] = (m_minimize) ? -c : c;

      // reset work parameter
      m_workParams[i % dims] = m_parameters[i % dims];

      // check if result is better than the last one
      if (i)
      {
        if (currvalue < m_workResults[i])
        {
          currvalue = m_workResults[i];
          bestpos = i;
        }
      }
      else
        currvalue = m_workResults[i];
      m_numEvaluations++;
    }

    // if best result is better, set it
    if (currvalue > lastvalue)
    {

      // set this direction permanently
      m_parameters[bestpos % dims] +=
        m_stepsizes[bestpos % dims] * m_scale * (bestpos <
            dims ? 1.0 : -1.0);
      m_workParams[bestpos % dims] = m_parameters[bestpos % dims];
      m_costFunctionValue = (m_minimize) ? -currvalue : currvalue;
      better = true;
      downscaled = false;
    }

    // otherwise scale the step size down
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


