
#include <array2d.h>
#include <SimilarityMeasures.h>
#include <OptimizerPowellBrent.h>
#include <OptimizerGradientDescent.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
  
/**
 *   @example OptimizerExample.cpp
 *    This is an example of how to use the Optimizer Classes
 */ 

class TestOptimizer : public OptimizerGradientDescent
{
  public:
    virtual double costFunction (double *tmpParams)
    {
      int x = (int) tmpParams[0];
      if (x < 0)
        x = 0;
      if (x >= width)
        x = width - 1;
      int y = (int) tmpParams[1];
      if (y < 0)
        y = 0;
      if (y >= height)
        y = height - 1;
      int z = (int) tmpParams[2];
      image1->fill (0);
      image1->set (y - z, y + z + 127, x - z, x + z + 127, 60);
      sm.CalculateMeasures ();
      return sm.GetNormalizedCrossCorr ();
    }
    virtual bool init ()
    {
      m_maxIterations = 10000;
      m_maxEvaluations = 1000000;
      m_functionTolerance = 0.000001;
      m_paramsTolerance = 0.000001;
      m_learningRate = 1000;
      width = 256;
      height = 256;
      image1 = new hs::array2d < char >(height, width);
      image2 = new hs::array2d < char >(height, width);
      image2->fill (0);
      image2->set (64, 191, 64, 191, 60);
      sm.SetFirstImage (*image1, 0, 8);
      sm.SetSecondImage (*image2, 0, 8);
      sm.SetMutualInformation (false);
      sm.SetNormalizedCrossCorr (true);
      sm.InitInternalBuffer (8);
      return true;
    }
    virtual void verbose ()
    {
      printf ("[ ");
      for (int i = 0; i < m_numberOfParams; i++)
        printf ("%f ", m_parameters[i]);
      printf ("] -> %f\n", m_costFunctionValue);
    } private:int width;
    int height;
    hs::array2d < char >*image1;
    hs::array2d < char >*image2;
    hs::reg2d3d::SimilarityMeasures < char >sm;
};

class TestOptimizerSimple : public OptimizerPowellBrent
{
  public:
    virtual double costFunction (double *tmpParams)
    {
      double val = 0, tmp;
      for (int i = 0; i < m_numberOfParams; i++)
      {

        //                      val += (tmpParams[i]-i)*(tmpParams[i]-i);
        if (i)
          tmp = tmpParams[i] - tmpParams[i - 1] - i;
        val += i ? tmp * tmp : tmpParams[i] * tmpParams[i];
      }
      return val;
    }
    virtual void verbose ()
    {
      printf ("[ ");
      for (int i = 0; i < m_numberOfParams; i++)
        printf ("%f ", m_parameters[i]);
      printf ("] -> %f\n", m_costFunctionValue);
    } };

void
OptimizerExample ()  {
  using namespace std;

  const int dims = 6;
  double params[dims] = { -2, -2, 4, 2, 1, 2 };
  double steps[dims] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  TestOptimizerSimple to;
  to.setInitialParameters (params, steps, dims);

  //      to.setLearningRate(0.8);

  /*	FILE *file = fopen("c:\\plot.txt", "w");
        for (int i = 48; i <= 80 ; i+=2) {
        for (int j = -16; j <= 16; j+=2) {
        params[0] = i; params[1] = 1; params[2] = j;
        double result = to.costFunction(params);
        fprintf(file, "%i %i %f\n", i, j, result);
        cout << result << endl;
        }
        fprintf(file, "\n");
        }
        printf("diagram written.\n");
        fclose(file); */ 
  to.optimize ();
  double *x = to.getParameters ();
  printf ("Result: [%f, %f, %f], f()=%f\n", x[0], x[1], x[2],
      to.getCostFunctionValue ());
  printf ("%i iterations, %i evaluations\n",
      to.getFinalNumberOfIterations (),
      to.getFinalNumberOfEvaluations ());
} 
