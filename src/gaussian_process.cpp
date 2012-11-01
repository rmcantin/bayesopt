/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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

GaussianProcess::GaussianProcess( double noise ):
  NonParametricProcess(noise)
{}  // Constructor


GaussianProcess::~GaussianProcess()
{} // Default destructor



double GaussianProcess::negativeLogLikelihood()
{
  const matrixd K = computeCorrMatrix();
  const size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  vectord alpha(mGPY);
  cholesky_solve(L,alpha,ublas::lower());
  double loglik = .5*ublas::inner_prod(mGPY,alpha) + trace(L) 
    + n*0.91893853320467; 
  // 0.9183... = log(2*pi)/2

  return loglik;
}


int GaussianProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  const double kn = (*mKernel)(query, query);
  const vectord kStar = computeCrossCorrelation(query);
  yPred = ublas::inner_prod(kStar,mAlphaV);

  vectord vd(kStar);
  ublas::inplace_solve(mL,vd,ublas::lower_tag());
  sPred = sqrt(kn - ublas::inner_prod(vd,vd));
  
  return 1;
}

ProbabilityDistribution* GaussianProcess::prediction(const vectord &query)
{
  double y, s;
  prediction(query,y,s);
  return new GaussianDistribution(y,s);
}


int GaussianProcess::precomputePrediction()
{
  const size_t n = mGPY.size();
  
  mAlphaV.resize(n,false);
  mAlphaV = mGPY;
  cholesky_solve(mL,mAlphaV,ublas::lower());

  return 1; 
}
	
