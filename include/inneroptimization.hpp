#ifndef __INNEROPTIMIZATION_HPP__
#define __INNEROPTIMIZATION_HPP__

#include "specialtypes.hpp"

// We plan to add more in the future since nlopt actually support many of them
enum innerOptAlgorithms {
  direct, lbfgs, combined
};


class InnerOptimization
{
public:
  InnerOptimization()
  { alg = direct; };
  virtual ~InnerOptimization(){};

  virtual double innerEvaluate(const vectord& query, 
			       vectord& grad)
  {return 0.0;}

  void setAlgorithm(innerOptAlgorithms newAlg)
  { alg = newAlg; }


protected:

  int innerOptimize(vectord &Xnext);
  int innerOptimize( double* x, int n, void* objPointer);	

  innerOptAlgorithms alg;
};

#endif
