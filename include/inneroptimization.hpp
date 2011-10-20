#ifndef __INNEROPTIMIZATION_HPP__
#define __INNEROPTIMIZATION_HPP__

#include "specialtypes.hpp"

class InnerOptimization
{
public:
  InnerOptimization(){};
  virtual ~InnerOptimization(){};

  virtual double innerEvaluate(const vectord& query)
  {return 0.0;}

protected:

  int innerOptimize(vectord &Xnext);
  int innerOptimize( double* x, int n, void* objPointer);	
};

#endif
