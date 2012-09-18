#include "bayesoptwpr.h"
#include "bayesoptcont.hpp"      


class CBayesOptContinuous: public BayesOptContinuous 
{
 public:

  CBayesOptContinuous( bopt_params params):
    BayesOptContinuous(params)
  {}; 


  double evaluateSample( const vectord &Xi ) 
  {
    int n = static_cast<int>(Xi.size());

    return  mF(n,&Xi[0],NULL,mOtherData);
  };

  void set_eval_funct(eval_func f)
  {  mF = f; }


  void save_other_data(void* other_data)
  {  mOtherData = other_data; }
 
protected:
  void* mOtherData;
  eval_func mF;
};



int bayes_optimization(int nDim, eval_func f, void* f_data,
		       const double *lb, const double *ub, /**< bounds */
		       double *x, /**< in: initial guess, out: minimizer */
		       double *minf, /**< out: minimum */
		       bopt_params parameters)
{
  vectord result(nDim);

  vectord lowerBound(nDim); 
  vectord upperBound(nDim); 

  std::copy(lb, lb+nDim, lowerBound.begin());
  std::copy(ub, ub+nDim, upperBound.begin());

  CBayesOptContinuous optimizer(parameters);

  optimizer.set_eval_funct(f);
  optimizer.save_other_data(f_data);
  optimizer.setBoundingBox(lowerBound,upperBound);
  optimizer.optimize(result);

  std::copy(result.begin(), result.end(), x);

  return 1; /* everything ok*/
}
