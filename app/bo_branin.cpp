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

#define _USE_MATH_DEFINES
#include <cmath>
#include <valarray>
#include "bayesoptcont.hpp"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

class ExampleBranin: public bayesopt::ContinuousModel
{
public:
  ExampleBranin(size_t dim,bopt_params par):
    ContinuousModel(dim,par) {}

  double evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }

    double x = xin(0) * 15 - 5;
    double y = xin(1) * 15;

    return sqr(y-(5.1/(4*sqr(M_PI)))*sqr(x)+5*x/M_PI-6)+10*(1-1/(8*M_PI))*cos(x)+10;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(2);  
    sv(0) = 0.1239; sv(1) = 0.8183;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.5428; sv(1) = 0.1517;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.9617; sv(1) = 0.1650;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

};

int main(int nargs, char *args[])
{
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = 200;
  par.n_init_samples = 50;
  par.theta[0] = 1.0;
  par.n_theta = 1;
  par.c_name = C_GP_HEDGE;
  
  ExampleBranin branin(2,par);
  vectord result(2);

  branin.optimize(result);
  std::cout << "Result:" << result << std::endl;
  branin.printOptimal();

  return 0;
}
