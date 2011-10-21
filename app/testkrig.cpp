#include <ctime>
#include "krigwpr.h"
#include "krigging.hpp"


class TestEGO: public SKO 
{
 public:

  double evaluateSample( const vectord &Xi ) 
  {
    double f = 10.;
    
    for (size_t i = 0; i < Xi.size(); ++i)
      {
	f += (Xi(i) - .53f) * (Xi(i) - .53f);
      }
    return f;
  };


  bool checkReachability( const vectord &query )
  {
    /*
      for (size_t i = 0; i < query.size(); ++i)
      {
	if ((query(i) > .6f) && (query(i) < .9f))
	  return false;
	  }*/
 
   return true;
  };
 
};

double testFunction(unsigned int n, double *x,
		    double *gradient, /* NULL if not needed */
		    void *func_data)
{/*
  double f = 10.;
    
  for (unsigned int i = 0; i < n; ++i)
    {
      f += (x[i] - .53f) * (x[i] - .53f);
      }*/

  double f = (x[0]-0.5)*(x[0]-0.5) + sin(20*x[0])*0.2;

  return f;
}



int main(int nargs, char *args[])
{    
  clock_t start, end;
  double diff,diff2;
  int n = 1;
  vectord result(n);
  randEngine rEng(100u);

  int nIterations = 600;

  
  std::cout << "Running C++ interface" << std::endl;
  TestEGO gp;
  start = clock();
  gp.optimize(result,rEng);
  end = clock();
  diff = (double)(end-start) / (double)CLOCKS_PER_SEC;


  std::cout << "Running C inferface" << std::endl;
  
  double u[128], x[128], l[128];
  double fmin[128];


  for (int i = 0; i < n; ++i) {
    l[i] = 0.;
    u[i] = 1.;
  }

  gp_params par;

  par.theta = KERNEL_THETA;
  par.p = KERNEL_P;
  par.alpha = PRIOR_ALPHA;
  par.beta = PRIOR_BETA;
  par.delta = PRIOR_DELTA_SQ;
  par.noise = DEF_REGULARIZER;
  
  start = clock();

  krigging_optimization(n,&testFunction,NULL,l,u,x,fmin,nIterations,par,1,0);

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

