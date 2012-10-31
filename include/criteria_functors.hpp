/** -*- c++ -*- \file criteria_functors.hpp \brief Criteria functions */
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

#ifndef  _CRITERIA_FUNCTORS_HPP_
#define  _CRITERIA_FUNCTORS_HPP_

#include "parameters.h"
#include "nonparametricprocess.hpp"

///\addtogroup CriteriaFunctions
//@{

/**
 * \brief Abstract interface for criteria functors.
 */
class Criteria
{
public:
  Criteria(NonParametricProcess *proc):  mProc(proc) {};
  virtual ~Criteria(){};
  virtual double operator()( const vectord &x) = 0;
  virtual void resetAnneal() {};  //dummy function

  static Criteria* create(criterium_name name,
			  NonParametricProcess* proc);
protected:
  NonParametricProcess *mProc;
};


/**
 * \brief Expected improvement criterion by Mockus.
 */
class ExpectedImprovement: public Criteria
{
public:
  ExpectedImprovement(NonParametricProcess *proc): Criteria(proc), mExp(1) {};
  virtual ~ExpectedImprovement(){};
  inline void setExponent(size_t exp) {mExp = exp;};
  double operator()(const vectord &x)
  { return mProc->negativeExpectedImprovement(x,mExp); };
private:
  size_t mExp;
};


/**
 * \brief Lower (upper) confidence bound criterion.
 */
class LowerConfidenceBound: public Criteria
{
public:
  LowerConfidenceBound(NonParametricProcess *proc):Criteria(proc),mBeta(1) {};
  virtual ~LowerConfidenceBound(){};
  inline void setBeta(double beta) { mBeta = beta; };
  double operator()( const vectord &x)
  { return mProc->lowerConfidenceBound(x,mBeta); };
private:
  double mBeta;
};


/**
 * \brief Probability of improvement criterion based on (Kushner).
 */
class ProbabilityOfImprovement: public Criteria
{
public:
  ProbabilityOfImprovement(NonParametricProcess *proc): 
    Criteria(proc), mEpsilon(0.01) {};
  virtual ~ProbabilityOfImprovement(){};
  inline void setEpsilon(double eps) {mEpsilon = eps;};
  double operator()( const vectord &x)
  { return mProc->negativeProbabilityOfImprovement(x,mEpsilon); };
private:
  double mEpsilon;
};

/**
 * \brief Greedy A-Optimality criterion. 
 * Used for learning the function, not to minimize.
 */
class GreedyAOptimality: public Criteria
{
public:
  GreedyAOptimality(NonParametricProcess *proc): Criteria(proc){};
  virtual ~GreedyAOptimality(){};
  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(x,yPred,sPred);
    return sPred;
  };
};

/**
 * \brief Expected return criterion.
 */
class ExpectedReturn: public Criteria
{
public:
  ExpectedReturn(NonParametricProcess *proc): Criteria(proc){};
  virtual ~ExpectedReturn(){};
  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(x,yPred,sPred);
    return yPred;
  };
};

/**
 * \brief Optimistic sampling, also called, Thomson sampling.
 */
class OptimisticSampling: public Criteria
{
public:
  OptimisticSampling(NonParametricProcess *proc): 
    Criteria(proc), mtRandom(100u) {};
  virtual ~OptimisticSampling(){};
  double operator()( const vectord &x)
  {
    double yPred, sPred, yStar = mProc->sample_query(x,mtRandom);
    mProc->prediction(x,yPred,sPred);
    return std::min(yPred,yStar);
  };
private:
  randEngine mtRandom;
};


/**
 * \brief Expected improvement criterion using Schonlau annealing.
 */
class AnnealedExpectedImprovement: public Criteria
{
public:
  AnnealedExpectedImprovement(NonParametricProcess *proc):
    Criteria(proc){ resetAnneal(); };
  virtual ~AnnealedExpectedImprovement(){};
  inline void setExponent(size_t exp) {mExp = exp;};
  void resetAnneal() { nCalls = 0; mExp = 10;};
  double operator()( const vectord &x)
  {
    ++nCalls;
    if (nCalls % 10)
      mExp = ceil(mExp/2.0);

    return mProc->negativeExpectedImprovement(x,mExp);
  };

private:
  size_t mExp;
  unsigned int nCalls;
};


/**
 * \brief Lower (upper) confidence bound using Srinivas annealing
 */
class AnnealedLowerConfindenceBound: public Criteria
{
public:
  AnnealedLowerConfindenceBound(NonParametricProcess *proc):
    Criteria(proc){ resetAnneal(); };
  virtual ~AnnealedLowerConfindenceBound(){};
  inline void setBetaCoef(double betac) { mCoef = betac; };
  void resetAnneal() { nCalls = 0; mCoef = 5.0;};
  double operator()( const vectord &x)
  {
    ++nCalls;
    size_t nDims = x.size();
    
    double beta = sqrt(2*log(static_cast<double>(nCalls*nCalls))*(nDims+1) 
		       + log(static_cast<double>(nDims))*nDims*mCoef);

    return mProc->lowerConfidenceBound(x,beta);
  };
private:
  double mCoef;
  unsigned int nCalls;
};

//@}

#endif
