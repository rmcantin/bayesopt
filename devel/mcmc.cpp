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
#include "mcmc.hpp"

namespace bayesopt
{
  MCMCSampler::MCMCSampler(RBOptimizable* rbo, size_t dim)
  {
    obj = new RBOptimizableWrapper(rbo);

    mAlg = SLICE_MCMC;
    mDims = dim;
    nBurnOut = 500;
    nSamples = 100;
    mStepOut = true;
    mSigma = svectord(dim,6);
  };

  MCMCSampler::~MCMCSampler()
  {
    if (obj != NULL) delete obj;
  };


  void MCMC::sliceSample(vectord &x)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    size_t n = x.size();

    std::vector<int> perms = utils::return_index_vector(n);
    utils::randomPerms(perms, mtRandom);

    for (size_t i = 0; i<n; ++i)
      {
	size_t ind = perms[i];
	double sigma = mSigma(ind);

	double y_max = obj->evaluate(x);
	double y = sample()*y_max;

	// Step out
	double x_cur = x(ind);
	double r = sample();
	double xl = x_cur - r * sigma;
	double xl = x_cur + (1-r)*sigma;

	if (mStepOut)
	  {
	    x(ind) = xl;
	    while (obj->evaluate(x) > y) x(ind) -= sigma;
	    xl = x(ind);

	    x(ind) = xr;
	    while (obj->evaluate(x) > y) x(ind) += sigma;
	    xr = x(ind);
	  }

	//Shrink
	bool on_slice = false;
	while (!on_slice)
	  {
	    x(ind) = (xr-xl) * sample() + xl;
	    if (obj->evaluate(x) < y)
	      {
		if      (x(ind) > x_cur)  xr = x(ind);
		else if (x(ind) < x_cur)  xl = x(ind);
		else throw std::runtime_error("Error in MCMC. Slice colapsed.");
	      }
	    else
	      {
		on_slice = true;
	      }
	  }
      }
  }

  //TODO: Include new algorithms when we add them.
  void MCMC::run(const vectord &initX)
  {
    vectord x = initX;
    if (nBurnOut<0) burnOut(x);
  
    mParticles.clear();
    for(size_t i=0; i<nSamples; ++i)  
      {
	sliceSample(x);
	mParticles.push_back(x);
      }
  }


}// namespace bayesopt
