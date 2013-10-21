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

#include "log.hpp"
#include "fullbayesprocess.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas; 

  FullBayesProcess::FullBayesProcess(size_t dim, bopt_params params):
    NonParametricProcess(dim,params),mGeneralParams(params),mWeights(N_PROC)
  {
    initializeKernelParameters();
  };

  FullBayesProcess::~FullBayesProcess(){};

  ProbabilityDistribution* FullBayesProcess::prediction(const vectord &query);
  {
    //Sum of Gaussians?
  };

  int FullBayesProcess::initializeKernelParameters()
  {
    double w = 1.0/static_cast<double>(N_PROC);
    mWeights = svectord(n,w);

    //All the inner processes share the same parameters except the
    //kernel paramenters and the learning type.
    bopt_params newParams = mGeneralParams;
    newParams.learning_type = L_FIXED;
    size_t nhp = mGeneralParams.kernel.n_hp;
    matrixd kTheta(nhp,N_PROC);
    randEngine reng(200u);
    lhs(kTheta,reng);

    for (size_t ii = 0; ii < N_PROC; ++ii)
      { 
	vectord th = column(kTheta,ii);
	std::copy(th.begin(),th.end(),newParams.kernel.hp_mean);
	mVProc.push_back(NonParametricProcess::create(dim_,newParams));
      }
    

    return 0;
  }

  int FullBayesProcess::updateKernelParameters()
  {
    double sum = 0.0;
    for (size_t ii = 0; ii < N_PROC; ++ii)
      { 
	double lik = mVProc.evaluateKernelParams();
	mWeights(ii) *= lik;
	sum += mWeights(ii);
      }
    mWeights /= sum;  //Normalization
  };
