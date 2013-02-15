#include "bayesoptwpr.h"
#include "bayesoptcont.hpp"      
#include "bayesoptdisc.hpp"      

/**
 * \brief Version of BayesOptContinuous for the C wrapper
 */
class CBayesOptContinuous: public BayesOptContinuous 
{
 public:

  CBayesOptContinuous(size_t dim, bopt_params params):
    BayesOptContinuous(dim,params)
  {}; 

  virtual ~CBayesOptContinuous(){};

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

/**
 * \brief Version of BayesOptDiscrete for the C wrapper
 */
class CBayesOptDiscrete: public BayesOptDiscrete
{
 public:

  CBayesOptDiscrete(const vecOfvec &validX, bopt_params params):
    BayesOptDiscrete(validX, params)
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
		       const double *lb, const double *ub,
		       double *x, double *minf, bopt_params parameters)
{
  vectord result(nDim);

  vectord lowerBound(nDim); 
  vectord upperBound(nDim); 

  std::copy(lb, lb+nDim, lowerBound.begin());
  std::copy(ub, ub+nDim, upperBound.begin());

  CBayesOptContinuous optimizer(nDim, parameters);

  optimizer.set_eval_funct(f);
  optimizer.save_other_data(f_data);
  optimizer.setBoundingBox(lowerBound,upperBound);
  optimizer.optimize(result);

  std::copy(result.begin(), result.end(), x);

  return 1; /* everything ok*/
};

int bayes_optimization_disc(int nDim, eval_func f, void* f_data,
			    double *valid_x, size_t n_points,
			    double *x, double *minf, bopt_params parameters)
{
  vectord result(nDim);
  vectord input(nDim);
  vecOfvec xSet;

  for(size_t i = 0; i<n_points;++i)
    {
      for(int j = 0; j<nDim; ++j)
	{
	 input(j) = valid_x[i*nDim+j]; 
	}
      xSet.push_back(input);
    }

  if(parameters.n_init_samples > n_points)
    {
      parameters.n_init_samples = n_points;
      parameters.n_iterations = 0;
    }

  CBayesOptDiscrete optimizer(xSet,parameters);

  optimizer.set_eval_funct(f);
  optimizer.save_other_data(f_data);
  optimizer.optimize(result);

  std::copy(result.begin(), result.end(), x);

  return 1; /* everything ok*/
}
