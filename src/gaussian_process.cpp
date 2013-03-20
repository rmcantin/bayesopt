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
#include "gaussian_process.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"
#include "gauss_distribution.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas;

  GaussianProcess::GaussianProcess(size_t dim, double noise):
    NonParametricProcess(dim, noise)
  {
    d_ = new GaussianDistribution();
  }  // Constructor


  GaussianProcess::~GaussianProcess()
  {
    delete d_;
  } // Default destructor



  double GaussianProcess::negativeLogLikelihood()
  {
    const matrixd K = computeCorrMatrix();
    const size_t n = K.size1();
    matrixd L(n,n);
    utils::cholesky_decompose(K,L);

    vectord alpha(mGPY-mMeanV);
    inplace_solve(L,alpha,ublas::lower_tag());
    double loglik = .5*ublas::inner_prod(alpha,alpha) + log_trace(L);
    return loglik;
  }

  ProbabilityDistribution* GaussianProcess::prediction(const vectord &query)
  {
    const double kn = (*mKernel)(query, query);
    const vectord kStar = computeCrossCorrelation(query);
    double yPred = ublas::inner_prod(kStar,mAlphaV);

    vectord vd(kStar);
    ublas::inplace_solve(mL,vd,ublas::lower_tag());
    double sPred = sqrt(kn - ublas::inner_prod(vd,vd));

    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }


  int GaussianProcess::precomputePrediction()
  {
    const size_t n = mGPY.size();
  
    mAlphaV.resize(n,false);
    mAlphaV = mGPY-mMeanV;
    utils::cholesky_solve(mL,mAlphaV,ublas::lower());

    return 1; 
  }
	
} //namespace bayesopt
