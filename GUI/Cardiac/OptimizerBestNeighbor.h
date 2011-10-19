/**
*  @file	BestNeighbor.h
*  @brief	Simple Best Neighbor Optimizer
*  @author	Hari Sundar
*  @date	6/2/03
*
**/  
  
#ifndef BN_OPTIMIZER_H
#define BN_OPTIMIZER_H
  
#include <Optimizer.h>

class OptimizerBestNeighbor:public Optimizer 
{
  public:
    virtual bool init ();
    virtual bool optimize ();
};

#endif	/*  */

