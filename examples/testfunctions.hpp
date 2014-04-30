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
#include <algorithm>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt.hpp"

class ExampleOneD: public bayesopt::ContinuousModel
{
public:
  ExampleOneD(size_t dim, bopt_params par):
    ContinuousModel(dim,par) {}

  double evaluateSample(const vectord& xin)
  {
    if (xin.size() != 1)
      {
	std::cout << "WARNING: This only works for 1D inputs." << std::endl
		  << "WARNING: Using only first component." << std::endl;
      }

    double x = xin(0);
    return (x-0.3)*(x-0.3) + sin(20*x)*0.2;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  void printOptimal()
  {
    std::cout << "Optimal:" << 0.23719 << std::endl;
  }

};


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

    return sqr(y-(5.1/(4*sqr(M_PI)))*sqr(x)
	       +5*x/M_PI-6)+10*(1-1/(8*M_PI))*cos(x)+10;
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

class ExampleHartmann6: public bayesopt::ContinuousModel
{
public:
  ExampleHartmann6(size_t dim,bopt_params par):
    ContinuousModel(dim,par), mA(4,6), mC(4), mP(4,6)
  {
    mA <<= 10.0,   3.0, 17.0,   3.5,  1.7,  8.0,
      0.05, 10.0, 17.0,   0.1,  8.0, 14.0,
      3.0,   3.5,  1.7,  10.0, 17.0,  8.0,
      17.0,   8.0,  0.05, 10.0,  0.1, 14.0;
    
    mC <<= 1.0, 1.2, 3.0, 3.2;

    mP <<= 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886,
      0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991,
      0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650,
      0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381;
  }

  double evaluateSample( const vectord& xin)
  {
    double y = 0.0;
    for(size_t i=0; i<4; ++i)
      {
	double sum = 0.0;
	for(size_t j=0;j<6; ++j)
	  {
	    double val = xin(j)-mP(i,j);
	    sum -= mA(i,j)*val*val;
	  }
	y -= mC(i)*std::exp(sum);
      }
    return y;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(6);
    sv <<= 0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573;

    std::cout << "Solution: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

private:
  matrixd mA;
  vectord mC;
  matrixd mP;
};
