/**  \file mcmcm.hpp \brief Markov Chain Monte Carlo algorithms */
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


#ifndef  _MCMC_HPP_
#define  _MCMC_HPP_

#include <boost/ptr_container/ptr_vector.hpp>
#include "randgen.hpp"
#include "optimizable.hpp"
#include "posteriormodel.hpp"

namespace bayesopt {

  // We plan to add more in the future 
  typedef enum {
    SLICE_MCMC           ///< Slice sampling
  } McmcAlgorithms;


  class MCMCSampler
  {
  public:
    MCMCSampler(RBOptimizable* rbo, size_t dim, randEngine& eng);
    virtual ~MCMCSampler();

    /** Sets the optimization algorithm  */
    void setAlgorithm(McmcAlgorithms newAlg);

    void setNParticles(size_t nParticles);

    void setNBurnOut(size_t nParticles);

    /** Compute the inner optimization algorithm 
     * @param Xnext input: initial point of the Markov Chain, 
     *              output: last point of the MC
     */
    void run(vectord &Xnext);

    vectord getParticle(size_t i);

    void printParticles();

  private:
    void randomJump(vectord &x);
    void burnOut(vectord &x);
    void sliceSample(vectord &x);

    RBOptimizableWrapper* obj;

    McmcAlgorithms mAlg;
    size_t mDims;
    size_t nBurnOut;
    size_t nSamples;
    bool mStepOut;

    vectord mSigma;
    vecOfvec mParticles;
    randEngine mtRandom;

  private: //Forbidden
    MCMCSampler();
    MCMCSampler(MCMCSampler& copy);
  };

  inline void MCMCSampler::randomJump(vectord &x)
  {
    randNFloat sample( mtRandom, normalDist(0,1) );
    FILE_LOG(logERROR) << "Doing random jump.";
    for(vectord::iterator it = x.begin(); it != x.end(); ++it)
      {
	*it = sample()*6;
      }
    FILE_LOG(logERROR) << "Likelihood." << x << " | " << obj->evaluate(x);
  }

  //TODO: Include new algorithms when we add them.
  inline void MCMCSampler::burnOut(vectord &x)
  {
    for(size_t i=0; i<nBurnOut; ++i)  
      {
	try
	  {
	    sliceSample(x);
	  }
	catch(std::runtime_error& e)
	  {
	    FILE_LOG(logERROR) << e.what();
	    randomJump(x);
	  }
      }
  }

  inline void MCMCSampler::setNParticles(size_t nParticles)
  { nSamples = nParticles; };

  inline void MCMCSampler::setNBurnOut(size_t nParticles)
  { nBurnOut = nParticles; };

  inline void MCMCSampler::setAlgorithm(McmcAlgorithms newAlg)
  { mAlg = newAlg; };

  inline vectord MCMCSampler::getParticle(size_t i)
  { return mParticles[i]; };

  inline void MCMCSampler::printParticles()
  {
    for(size_t i=0; i<mParticles.size(); ++i)
      std::cout << i << "->" << mParticles[i] << std::endl;
  }


  class MCMCModel: public PosteriorModel
  {
  public:

    typedef boost::ptr_vector<NonParametricProcess>  GPVect;
    typedef boost::ptr_vector<Criteria>  CritVect;

    /** 
     * Constructor
     * @param params set of parameters (see parameters.h)
     */
    MCMCModel(size_t dim, bopt_params params, randEngine& eng);

    /** 
     * Default destructor
     */
    virtual ~MCMCModel();

    void updateHyperParameters();
    void fitSurrogateModel();
    void updateSurrogateModel();
    double evaluateCriteria(const vectord& query);

    bool criteriaRequiresComparison();
    void setFirstCriterium();
    bool setNextCriterium(const vectord& prevResult);
    std::string getBestCriteria(vectord& best);

    Criteria* getCriteria();
    NonParametricProcess* getSurrogateModel();
   
  private:

    MCMCModel();

    void setSurrogateModel(randEngine& eng);    
    void setCriteria(randEngine& eng);

  private:  // Members
    size_t nParticles;
    GPVect mGP;                ///< Pointer to surrogate model
    CritVect mCrit;                    ///< Metacriteria model

    boost::scoped_ptr<MCMCSampler> kSampler;
  };

  /**@}*/

  inline void MCMCModel::fitSurrogateModel()
  { 
    for(GPVect::iterator it=mGP.begin(); it != mGP.end(); ++it)
      it->fitSurrogateModel(); 
  };

  inline void MCMCModel::updateSurrogateModel()
  {     
    for(GPVect::iterator it=mGP.begin(); it != mGP.end(); ++it)
      it->updateSurrogateModel(); 
  };

  inline double MCMCModel::evaluateCriteria(const vectord& query)
  { 
    double sum = 0.0;
    for(CritVect::iterator it=mCrit.begin(); it != mCrit.end(); ++it)
      {
	sum += it->evaluate(query); 
      }
    return sum/static_cast<double>(nParticles);
  };

  inline bool MCMCModel::criteriaRequiresComparison()
  {return mCrit[0].requireComparison(); };
    
  inline void MCMCModel::setFirstCriterium()
  { 
    for(CritVect::iterator it=mCrit.begin(); it != mCrit.end(); ++it)
      {
	it->initialCriteria();
      }
  };

  // Although we change the criteria for all MCMC particles, we use
  // only the first element to compute de Hedge algorithm, because it
  // should be based on the average result, thus being common for all
  // the particles.
  inline bool MCMCModel::setNextCriterium(const vectord& prevResult)
  { 
    bool rotated;
    mCrit[0].pushResult(prevResult);
    for(CritVect::iterator it=mCrit.begin(); it != mCrit.end(); ++it)
      {
	rotated = it->rotateCriteria();
      }
    return rotated; 
  };

  inline std::string MCMCModel::getBestCriteria(vectord& best)
  { return mCrit[0].getBestCriteria(best); };


  inline  Criteria* MCMCModel::getCriteria()
  { return &mCrit[0]; };

  inline  NonParametricProcess* MCMCModel::getSurrogateModel()
  { return &mGP[0]; };



} //namespace bayesopt


#endif
