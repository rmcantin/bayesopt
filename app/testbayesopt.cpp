#include <ctime>
#include "bayesoptwpr.h"             // For the C-AP
#include "bayesoptcont.hpp"              // For the C++-API


/* Function to be used for C-API testing */
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

/* Class to be used for C++-API testing */
class TestEGO: public SKO 
{
 public:

  TestEGO(NonParametricProcess* gp):
    SKO(gp) {}

  double evaluateSample( const vectord &Xi ) 
  {
    double x[100];
    for (size_t i = 0; i < Xi.size(); ++i)
	x[i] = Xi(i);
    return testFunction(Xi.size(),x,NULL,NULL);
  };


  bool checkReachability( const vectord &query )
  { return true; };
 
};


int main(int nargs, char *args[])
{    
  int n = 1;                   // Number of dimensions
  int nIterations = 300;       // Number of iterations

  // Common configuration
  // See ctypes.h for the available options
  criterium_name c_name = c_ei;
  surrogate_name s_name = s_gaussianProcess;
  kernel_name k_name = k_seiso;
  gp_params par;

  par.theta = KERNEL_THETA;
  par.alpha = PRIOR_ALPHA;
  par.beta = PRIOR_BETA;
  par.delta = PRIOR_DELTA_SQ;
  par.noise = DEF_REGULARIZER;
  /*******************************************/

  clock_t start, end;
  double diff,diff2;

  std::cout << "Running C++ interface" << std::endl;
  NonParametricProcess* gp;

  // Configure C++ interface
  switch(s_name)
    {
    case s_gaussianProcess: 
      gp = new BasicGaussianProcess(par.theta,par.noise);
      break;
    case s_gaussianProcessHyperPriors: 
      gp = new GaussianProcess(par.theta,par.noise,
			       par.alpha,par.beta,
			       par.delta);
      break;
    case s_studentTProcess:
      gp = new StudentTProcess(par.theta,par.noise);
      break; 
    }
  gp->setKernel(k_name);

  TestEGO gp_opt(gp);
  vectord result(n);
  gp_opt.setCriteria(c_name);

  // Run C++ interface
  start = clock();
  gp_opt.optimize(result,nIterations);
  end = clock();
  diff = (double)(end-start) / (double)CLOCKS_PER_SEC;
  /*******************************************/


  std::cout << "Running C inferface" << std::endl;
  
  // Configure C interface
  double u[128], x[128], l[128], fmin[128];

  for (int i = 0; i < n; ++i) {
    l[i] = 0.;    u[i] = 1.;
  }

  // Run C interface
  start = clock();
  bayes_optimization(n,&testFunction,NULL,l,u,x,fmin,
		     nIterations,par,c_name,s_name,k_name);
  end = clock();
  diff2 = (double)(end-start) / (double)CLOCKS_PER_SEC;
  /*******************************************/


  // Results
  std::cout << "Final result C++: " << result << std::endl;
  std::cout << "Final result C: (";
  for (int i = 0; i < n; i++ )
    std::cout << x[i] << ", ";
  
  std::cout << ")" << std::endl;

  std::cout << "Elapsed time in C++: " << diff << " seconds" << std::endl;
  std::cout << "Elapsed time in C: " << diff2 << " seconds" << std::endl;
}

