#include "lhs.hpp"
#include "bayesoptdisc.hpp"
#include "bayesoptwpr.h"



/* Function to be used for C-API testing */
double testFunction(unsigned int n, const double *x,
		    double *gradient, /* NULL if not needed */
		    void *func_data)
{
  double f = 10.;
  for (unsigned int i = 0; i < n; ++i)
    {
      f += (x[i] - .53) * (x[i] - .53);
    }
  return f;
}

/* Class to be used for C++-API testing */
class TestDisc: public bayesopt::BayesOptDiscrete
{
 public:

  TestDisc(const vecOfvec & validSet, bopt_params param):
    bayesopt::BayesOptDiscrete(validSet,param) {}

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
  int n = 6;                   // Number of dimensions

  // Common configuration
  // See parameters.h for the available options
  // If we initialize the struct with the DEFAUL_PARAMS,
  // the we can optionally change only few of them 
  bopt_params par = initialize_parameters_to_default();

  par.theta[0] = KERNEL_THETA;
  par.n_theta = 1;
  par.alpha = PRIOR_ALPHA;
  par.beta = PRIOR_BETA;
  par.delta = PRIOR_DELTA_SQ;
  par.noise = DEFAULT_NOISE;
  par.c_name = C_EI;
  par.s_name = S_GAUSSIAN_PROCESS_INV_GAMMA_NORMAL;
  par.k_name = K_MATERN_ISO3;
  par.n_iterations = 20;       // Number of iterations
  par.n_init_samples = 20;
  /*******************************************/
  
  size_t nPoints = 50;

  randEngine mtRandom;
  matrixd xPoints(nPoints,n);
  vecOfvec xP;

  //size_t nPinArr = n*nPoints;
  double xPointsArray[1000];

  lhs(xPoints,mtRandom);

  for(size_t i = 0; i<nPoints; ++i)
    {
      vectord point = row(xPoints,i);  
      xP.push_back(point);
      for(size_t j=0; j<n; ++j)
	{
	  xPointsArray[i*n+j] = point(j);	  
	}
    }
    

  std::cout << "Running C++ interface" << std::endl;
  // Run C++ interface
  TestDisc gp_opt(xP,par);
  vectord result(n);
  gp_opt.optimize(result);

  std::cout << "Running C interface" << std::endl;
  // Run C interface

  double x[128], fmin[128];
  bayes_optimization_disc(n, &testFunction, NULL, xPointsArray, nPoints,
			  x, fmin, par);

  // Find the optimal value
  size_t min = 0;
  double minvalue = gp_opt.evaluateSample(row(xPoints,min));
  for(size_t i = 1; i<nPoints; ++i)
    {
      vectord point = row(xPoints,i);  
      if (gp_opt.evaluateSample(point) < minvalue)
	{
	  min = i;
	  minvalue = gp_opt.evaluateSample(row(xPoints,min));
	}
    }

  std::cout << "Final result C++: " << result << std::endl;
  std::cout << "Final result C: (";
  for (int i = 0; i < n; i++ )
    std::cout << x[i] << ", ";
  std::cout << ")" << std::endl;
  std::cout << "Optimal: " << row(xPoints,min) << std::endl;
}
