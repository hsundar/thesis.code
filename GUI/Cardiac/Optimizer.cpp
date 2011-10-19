/**
*  @file	Optimizer.cpp
*  @brief	Implements the optimizer class
*  @author	Hari Sundar
*  @date	5/29/03
*
**/  

#include "Optimizer.h"

#ifndef NULL
#define NULL 0
#endif

Optimizer::Optimizer () 
{
  m_numberOfParams = 0;
  m_parameters = NULL;
  m_stepsizes = NULL;
  m_workParams = NULL;
  m_workResults = NULL;
  m_minimize = false;	// default is maximizing
  m_isOptimizing = false;
}
Optimizer::~Optimizer () 
{
  if (m_parameters)
    delete[]m_parameters;
  if (m_workParams)
    delete[]m_workParams;
  if (m_workResults)
    delete[]m_workResults;
}
void Optimizer::setInitialParameters (double *parameters,
    double *stepsizes,
    int dimensions) 
{
  if (m_parameters)
    delete[]m_parameters;
  if (m_stepsizes)
    delete[]m_stepsizes;
  m_numberOfParams = dimensions;
  m_parameters = new double[dimensions];
  m_stepsizes = new double[dimensions];
  for (int i = 0; i < dimensions; i++)
  {
    m_parameters[i] = parameters[i];
    m_stepsizes[i] = stepsizes[i];
  } 
  init ();
};

