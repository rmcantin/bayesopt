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
#include "lhs.hpp"
#include "mcmc.hpp"

namespace bayesopt
{
  MCMCSampler::MCMCSampler(RBOptimizable* rbo, size_t dim, randEngine& eng):
    mtRandom(eng)
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


  void MCMCSampler::sliceSample(vectord &x)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    size_t n = x.size();

    std::vector<int> perms = utils::return_index_vector(0,n);

    utils::randomPerms(perms, mtRandom);

    for (size_t i = 0; i<n; ++i)
      {
	const size_t ind = perms[i];
	const double sigma = mSigma(ind);

	const double y_max = obj->evaluate(x);
	const double y = y_max-std::log(sample());  
	//y = y_max * sample(), but we are in negative log space

	if (y == 0.0) 
	  {
	    throw std::runtime_error("Error in MCMC: Initial point out of support region."); 
	  }

	// Step out
	const double x_cur = x(ind);
	const double r = sample();
	double xl = x_cur - r * sigma;
	double xr = x_cur + (1-r)*sigma;

	if (mStepOut)
	  {
	    x(ind) = xl;
	    while (obj->evaluate(x) < y) { x(ind) -= sigma; }
	    xl = x(ind);

	    x(ind) = xr;
	    while (obj->evaluate(x) < y) { x(ind) += sigma; }
	    xr = x(ind);
	  }

	//Shrink
	bool on_slice = false;
	while (!on_slice)
	  {
	    x(ind) = (xr-xl) * sample() + xl;
	    if (obj->evaluate(x) > y)
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
  void MCMCSampler::run(vectord &Xnext)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    if (nBurnOut>0) burnOut(Xnext);

    mParticles.clear();
    for(size_t i=0; i<nSamples; ++i)  
      {
	try
	  {
	    sliceSample(Xnext);
	  }
	catch(std::runtime_error& e)
	  {
	    FILE_LOG(logERROR) << e.what();
	    randomJump(Xnext);
	  }
	mParticles.push_back(Xnext);

      }
  }

  //////////////////////////////////////////////////////////////////////
  MCMCModel::MCMCModel(size_t dim, bopt_params parameters, 
		       randEngine& eng):
    PosteriorModel(dim,parameters,eng), nParticles(10)
  {
    //TODO: Take nParticles from parameters
    
    // Configure Surrogate and Criteria Functions
    setSurrogateModel(eng);
    setCriteria(eng);

    // Seting MCMC for kernel hyperparameters...
    // We use the first GP as the "walker" to get the particles. Then,
    // we will use a whole vector of GPs to avoid recomputing the
    // kernel matrices after every data point.
    size_t nhp = mGP[0].nHyperParameters();
    kSampler.reset(new MCMCSampler(&mGP[0],nhp,eng));

    kSampler->setNParticles(nParticles);
    kSampler->setNBurnOut(300);
  }

  MCMCModel::~MCMCModel()
  { } // Default destructor


  void MCMCModel::updateHyperParameters()
  {
    // We take the initial point as the last particle from previous update.
    size_t last = mGP.size()-1;
    vectord lastTheta = mGP[last].getHyperParameters();

    FILE_LOG(logDEBUG) << "Initial kernel parameters: " << lastTheta;
    kSampler->run(lastTheta);
    for(size_t i = 0; i<nParticles; ++i)
      {
	mGP[i].setHyperParameters(kSampler->getParticle(i));
      }
    FILE_LOG(logDEBUG) << "Final kernel parameters: " << lastTheta;	
  };


  void MCMCModel::setSurrogateModel(randEngine& eng)
  {
    for(size_t i = 0; i<nParticles; ++i)
      {
	mGP.push_back(NonParametricProcess::create(mDims,mParameters,
						   mData,mMean,eng));
      } 
  } // setSurrogateModel

  void MCMCModel::setCriteria(randEngine& eng)
  {
    CriteriaFactory mCFactory;

    for(size_t i = 0; i<nParticles; ++i)
      {
	mCrit.push_back(mCFactory.create(mParameters.crit_name,&mGP[i]));
	mCrit[i].setRandomEngine(eng);

	if (mCrit[i].nParameters() == mParameters.n_crit_params)
	  {
	    mCrit[i].setParameters(utils::array2vector(mParameters.crit_params,
						       mParameters.n_crit_params));
	  }
	else // If the number of paramerters is different, use default.
	  {
	    if (mParameters.n_crit_params != 0)
	      {
		FILE_LOG(logERROR) << "Expected " << mCrit[i].nParameters() 
				   << " parameters. Got " 
				   << mParameters.n_crit_params << " instead.";
	      }
	    FILE_LOG(logINFO) << "Usign default parameters for criteria.";
	  }
      }
  } // setCriteria




}// namespace bayesopt
