/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include <ctime>
#include "bayesoptwpr.h"                 // For the C API
#include "bayesoptcont.hpp"              // For the C++ API


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
class ExampleQuadratic: public bayesopt::ContinuousModel
{
 public:

  ExampleQuadratic(size_t dim,bopt_params param):
    ContinuousModel(dim,param) {}

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
  int n = 10;                   // Number of dimensions

  // Common configuration
  // See parameters.h for the available options
  // If we initialize the struct with the DEFAUL_PARAMS,
  // the we can optionally change only few of them 
  bopt_params par = initialize_parameters_to_default();

  par.kernel.theta[0] = KERNEL_THETA;
  par.kernel.theta[1] = KERNEL_THETA;
  par.kernel.s_theta[0] = 1;
  par.kernel.s_theta[1] = 1;
  par.kernel.n_theta = 2;
  par.mean.mu[0] = 1.0;
  par.mean.s_mu[0] = PRIOR_DELTA_SQ;
  par.mean.n_mu = 1;
  par.alpha = PRIOR_ALPHA;
  par.beta = PRIOR_BETA;
  par.noise = DEFAULT_NOISE;
  par.surr_name = S_STUDENT_T_PROCESS_JEFFREYS;
  par.kernel.name = "kSum(kSEISO,kConst)";
  par.mean.name = "mConst";
  par.n_iterations = 200;       // Number of iterations
  par.n_init_samples = 50;
  par.verbose_level = 2;
  /*******************************************/

  clock_t start, end;
  double diff,diff2;

  std::cout << "Running C++ interface" << std::endl;
  // Configure C++ interface

  ExampleQuadratic opt(n,par);
  vectord result(n);

  // Run C++ interface
  start = clock();
  opt.optimize(result);
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
  bayes_optimization(n,&testFunction,NULL,l,u,x,fmin,par);
  end = clock();
  diff2 = (double)(end-start) / (double)CLOCKS_PER_SEC;
  /*******************************************/


  // Results
  std::cout << "Final result C++: " << result << std::endl;
  std::cout << "Elapsed time in C++: " << diff << " seconds" << std::endl;

  std::cout << "Final result C: [" << n <<"](" << x[0];
  for (int i = 1; i < n; ++i )
    std::cout << "," << x[i];
  
  std::cout << ")" << std::endl;
  std::cout << "Elapsed time in C: " << diff2 << " seconds" << std::endl;

}

