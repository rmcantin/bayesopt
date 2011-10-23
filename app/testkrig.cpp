#include <ctime>
#include "krigwpr.h"
#include "krigging.hpp"

#include "gaussprocess.hpp"
#include "basicgaussprocess.hpp"
#include "studenttprocess.hpp"

double testFunction(unsigned int n, double *x,
		    double *gradient, /* NULL if not needed */
		    void *func_data)
{
  double f = 10.;
  for (unsigned int i = 0; i < n; ++i)
    {
      f += (x[i] - .53f) * (x[i] - .53f);
    }
  return f;
}



class TestEGO: public SKO 
{
 public:

  TestEGO( gp_params params, size_t nIter, 
	   NonParametricProcess* gp):
    SKO(params,nIter, gp) {}

  double evaluateSample( const vectord &Xi ) 
  {
    double x[100];
    for (size_t i = 0; i < Xi.size(); ++i)
	x[i] = Xi(i);
    return testFunction(Xi.size(),x,NULL,NULL);
  };


  bool checkReachability( const vectord &query )
  {
    return true;
  };
 
};


int main(int nargs, char *args[])
{    
  clock_t start, end;
  double diff,diff2;
  int n = 4;
  vectord result(n);

  int nIterations = 600;

  criterium_name c_name = c_ei;
  gp_params par;

  par.theta = KERNEL_THETA;
  par.alpha = PRIOR_ALPHA;
  par.beta = PRIOR_BETA;
  par.delta = PRIOR_DELTA_SQ;
  par.noise = DEF_REGULARIZER;
  
  std::cout << "Running C++ interface" << std::endl;
  NonParametricProcess* gp = new StudentTProcess(par.theta,par.noise);
  TestEGO gp_opt(par,nIterations,gp);

  start = clock();
  gp_opt.setCriteria(c_name);
  gp_opt.optimize(result);
  end = clock();
  diff = (double)(end-start) / (double)CLOCKS_PER_SEC;


  std::cout << "Running C inferface" << std::endl;
  
  double u[128], x[128], l[128];
  double fmin[128];


  for (int i = 0; i < n; ++i) {
    l[i] = 0.;
    u[i] = 1.;
  }
  start = clock();

  krigging_optimization(n,&testFunction,NULL,l,u,x,fmin,
			nIterations,par,c_name,s_studentTProcess);

  end = clock();
  diff2 = (double)(end-start) / (double)CLOCKS_PER_SEC;

  std::cout << "Final result C++: " << result << std::endl;
  std::cout << "Final result C: (";
  for (int i = 0; i < n; i++ )
    std::cout << x[i] << ", ";
  
  std::cout << ")" << std::endl;

  std::cout << "Elapsed time in C++: " << diff << " seconds" << std::endl;
  std::cout << "Elapsed time in C: " << diff2 << " seconds" << std::endl;
}

