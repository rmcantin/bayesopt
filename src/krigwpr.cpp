#include "krigwpr.h"
#include "krigging.hpp"



void copyCarrayInVector(const double* array, int n, vector<double>& vec)
{
  for(int i = 0; i<n;i++)
    vec(i) = array[i];
}


void copyVectorInArray(double* array, int n, const vector<double>& vec)
{
  for(int i = 0; i<n;i++)
    array[i] = vec(i);
}


class CKrigging: public Krigging 
{
 public:

  CKrigging( gp_params params,
	     size_t nIter, bool useCool = false): Krigging(params,nIter,useCool)
  {}; 


  double evaluateSample( const vector<double> &Xi ) 
  {
    int n = static_cast<int>(Xi.size());
    double x[n];

    copyVectorInArray(x,n,Xi);
     
    return mF(n,x,NULL,mOtherData);
  };

  void set_eval_funct(eval_func f)
  {
    mF = f;
  }


  void save_other_data(void* other_data)
  {
    mOtherData = other_data;
  }

 protected:

  void* mOtherData;
  eval_func mF;
};



int krigging_optimization(int nDim, eval_func f, void* f_data,
			  const double *lb, const double *ub, /* bounds */
			  double *x, /* in: initial guess, out: minimizer */
			  double *minf, /* out: minimum */
			  int maxeval, gp_params params,
			  int use_cooling_scheme)
{

  bool useCool = false;

  if (use_cooling_scheme)
    useCool = true;

  vector<double> result(nDim);
  randEngine rEng(100u);

  vector<double> lowerBound(nDim); 
  vector<double> upperBound(nDim); 

  copyCarrayInVector(lb,nDim,lowerBound);
  copyCarrayInVector(ub,nDim,upperBound);

  CKrigging optimizer(params, maxeval, useCool);

  optimizer.set_eval_funct(f);
  optimizer.save_other_data(f_data);
  optimizer.optimize(result,rEng);

  copyVectorInArray(x,nDim,result);

  return 1; /* everything ok*/
  
}
