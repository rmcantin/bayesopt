/**  \file criteria_atomic.hpp \brief Atomic (single) criterion functions */
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

#ifndef  _CRITERIA_ATOMIC_HPP_
#define  _CRITERIA_ATOMIC_HPP_

#include "criteria_functors.hpp"
#include "prob_distribution.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions
   * \brief Set of criterium functions to select the next point during
   * optimization/exploration
   */
  //@{

  /// Abstract class for atomic criteria (only one function)
  class AtomicCriteria: public Criteria
  {
  public:
    virtual ~AtomicCriteria(){};
    virtual int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      return 0;
    };
    bool requireComparison(){ return false; };
    virtual double operator()(const vectord &x)  = 0;
  };


  /// Expected improvement criterion by Mockus \cite Mockus78
  class ExpectedImprovement: public AtomicCriteria
  {
  public:
    virtual ~ExpectedImprovement(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mExp = 1;
      return 0;
    };

    int setParameters(const vectord &params)
    {
      mExp = static_cast<size_t>(params(0));
      return 0;
    };

    size_t nParameters() {return 1;};

    double operator()(const vectord &x)
    { 
      const double min = mProc->getValueAtMinimum();
      return mProc->prediction(x)->negativeExpectedImprovement(min,mExp); 
    };

    std::string name() {return "cEI";};

  private:
    size_t mExp;
  };

  /// Expected improvement criterion modification by Lizotte
  class BiasedExpectedImprovement: public AtomicCriteria
  {
  public:
    virtual ~BiasedExpectedImprovement(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mBias = 0.01;
      mExp = 1;
      return 0;
    };
    int setParameters(const vectord &params)
    {
      mExp = static_cast<size_t>(params(0));
      mBias = params(1);
      return 0;
    };

    size_t nParameters() {return 2;};

    double operator()(const vectord &x)
    { 
      const double sigma = mProc->getSignalVariance();
      const double min = mProc->getValueAtMinimum() - mBias/sigma;
      return mProc->prediction(x)->negativeExpectedImprovement(min,mExp); 
    };
    std::string name() {return "cBEI";};
  private:
    double mBias;
    size_t mExp;
  };


  /// Lower (upper) confidence bound criterion by [Cox and John, 1992].
  class LowerConfidenceBound: public AtomicCriteria
  {
  public:
    virtual ~LowerConfidenceBound(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mBeta = 1.0;
      return 0;
    };
    int setParameters(const vectord &params)
    {
      mBeta = params(0);
      return 0;
    };

    size_t nParameters() {return 1;};

    double operator()( const vectord &x)
    { 
      return mProc->prediction(x)->lowerConfidenceBound(mBeta); 
    };
    std::string name() {return "cLCB";};
  private:
    double mBeta;
  };


  /// Probability of improvement criterion based on (Kushner).
  class ProbabilityOfImprovement: public AtomicCriteria
  {
  public:
    virtual ~ProbabilityOfImprovement(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mEpsilon = 0.01;
      return 0;
    };
    int setParameters(const vectord &params)
    {
      mEpsilon = params(0);
      return 0;
    };

    size_t nParameters() {return 1;};

    inline void setEpsilon(double eps) { mEpsilon = eps; };
    double operator()( const vectord &x)
    { 
      const double min = mProc->getValueAtMinimum();
      return mProc->prediction(x)->negativeProbabilityOfImprovement(min,
								    mEpsilon); 
    };
    std::string name() {return "cPOI";};
  private:
    double mEpsilon;
  };

  /**
   * \brief Greedy A-Optimality criterion.  Used for learning the
   * function, not to minimize. Some authors name it I-optimality
   * because it minimizes the error on the prediction, not on the
   * parameters.
   */
  class GreedyAOptimality: public AtomicCriteria
  {
  public:
    virtual ~GreedyAOptimality(){};
    int setParameters(const vectord &params) { return 0; };
    size_t nParameters() {return 0;};
    double operator()( const vectord &x)
    { return mProc->prediction(x)->getStd(); };
    std::string name() {return "cAopt";};
  };


  /// Expected return criterion.
  class ExpectedReturn: public AtomicCriteria
  {
  public:
    virtual ~ExpectedReturn(){};
    int setParameters(const vectord &params) { return 0; };
    size_t nParameters() {return 0;};
    double operator()( const vectord &x)
    { return mProc->prediction(x)->getMean(); };
    std::string name() {return "cExpReturn";};
  };


  /**
   * \brief Optimistic sampling. A simple variation of Thompson sampling
   * that picks only samples that are better than the best outcome so
   * far.
   */
  class OptimisticSampling: public AtomicCriteria
  {
  public:
    OptimisticSampling(): mtRandom(100u) {};
    virtual ~OptimisticSampling(){};
    int setParameters(const vectord &params) { return 0; };
    size_t nParameters() {return 0;};
    double operator()( const vectord &x)
    {
      ProbabilityDistribution* d_ = mProc->prediction(x);
      const double yStar = d_->sample_query(mtRandom);
      const double yPred = d_->getMean();
      return (std::min)(yPred,yStar);
    };
    std::string name() {return "cOptimisticSampling";};
  private:
    randEngine mtRandom;
  };


  /**
   * \brief Thompson sampling. 
   * Picks a random realization of the surrogate model.
   */
  class ThompsonSampling: public AtomicCriteria
  {
  public:
    ThompsonSampling(): mtRandom(100u) {};
    virtual ~ThompsonSampling(){};
    int setParameters(const vectord &params) { return 0; };
    size_t nParameters() {return 0;};
    double operator()( const vectord &x)
    {
      ProbabilityDistribution* d_ = mProc->prediction(x);
      return d_->sample_query(mtRandom);
    };
    std::string name() {return "cThompsonSampling";};
  private:
    randEngine mtRandom;
  };



  /// Expected improvement criterion using Schonlau annealing. \cite Schonlau98
  class AnnealedExpectedImprovement: public AtomicCriteria
  {
  public:
    virtual ~AnnealedExpectedImprovement(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      reset();
      return 0;
    };

    int setParameters(const vectord &params)
    {
      mExp = static_cast<size_t>(params(0));
      return 0;
    };

    size_t nParameters() {return 1;};
    void reset() { nCalls = 0; mExp = 10;};
    double operator()( const vectord &x)
    {
      ++nCalls;
      if (nCalls % 10)
	mExp = ceil(mExp/2.0);

      ProbabilityDistribution* d_ = mProc->prediction(x);
      const double min = mProc->getValueAtMinimum();
      return d_->negativeExpectedImprovement(min,mExp); 
    };
    std::string name() {return "cEIa";};
  private:
    size_t mExp;
    unsigned int nCalls;
  };


  /// Lower (upper) confidence bound using Srinivas annealing \cite Srinivas10
  class AnnealedLowerConfindenceBound: public AtomicCriteria
  {
  public:
    virtual ~AnnealedLowerConfindenceBound(){};
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      reset();
      return 0;
    };

    int setParameters(const vectord &params)
    {
      mCoef = params(0);
      return 0;
    };
    size_t nParameters() {return 1;};
    void reset() { nCalls = 0; mCoef = 5.0;};
    double operator()( const vectord &x)
    {
      ++nCalls;
      size_t nDims = x.size();
    
      double beta = sqrt(2*log(static_cast<double>(nCalls*nCalls))*(nDims+1) 
			 + log(static_cast<double>(nDims))*nDims*mCoef);

      ProbabilityDistribution* d_ = mProc->prediction(x);
      return d_->lowerConfidenceBound(beta); 
    };
    std::string name() {return "cLCBa";};
  private:
    double mCoef;
    unsigned int nCalls;
  };

  /**
   * \brief Distance in input space. Can be combined with other
   * critera to trade off large changes in input space.
   */
  class InputDistance: public AtomicCriteria
  {
  public:
    int init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mW = 1;
      return 0;
    };
    virtual ~InputDistance(){};
    int setParameters(const vectord &params)
    {
      mW = params(0);
      return 0;
    };
    size_t nParameters() {return 1;};
 
    double operator()(const vectord &x)
    { 
      vectord x2(x.size());
      mProc->getLastSample(x2);
      return mW*norm_2(x-x2);
    };
    std::string name() {return "cDistance";};
  private:
    double mW;
  };


  //@}

} //namespace bayesopt


#endif
