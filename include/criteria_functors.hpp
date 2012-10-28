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
  Criteria(NonParametricProcess *proc):
    mProc(proc) {};

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
  ExpectedImprovement(NonParametricProcess *proc): 
    Criteria(proc), mExp(1) {};

  virtual ~ExpectedImprovement(){};
  inline void setExponent(size_t exp) {mExp = exp;};
  double operator()(const vectord &x)
  {
    double result = mProc->negativeExpectedImprovement(x);
    return result;
  };
private:
  size_t mExp;
};


/**
 * \brief Lower (upper) confidence bound criterion.
 */
class LowerConfidenceBound: public Criteria
{
public:
  LowerConfidenceBound(NonParametricProcess *proc): Criteria(proc) { mBeta = 1; };

  virtual ~LowerConfidenceBound(){};
  inline void setBeta(double beta) { mBeta = beta; };

  double operator()( const vectord &x)
  {
    return mProc->lowerConfidenceBound(x,mBeta);
  };
protected:
  double mBeta;
};


/**
 * \brief Probability of improvement criterion based on (Kushner).
 */
class ProbabilityOfImprovement: public Criteria
{
public:
  ProbabilityOfImprovement(NonParametricProcess *proc): Criteria(proc) {mEpsilon = 0.01;};
  virtual ~ProbabilityOfImprovement(){};
  inline void setEpsilon(double eps) {mEpsilon = eps;};

  double operator()( const vectord &x)
  {
    return mProc->negativeProbabilityOfImprovement(x,mEpsilon);
  };
protected:
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
    Criteria(proc), mtRandom(100u)
  {};

  virtual ~OptimisticSampling(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred, yStar = mProc->sample_query(x,mtRandom);
    mProc->prediction(x,yPred,sPred);
    return std::min(yPred,yStar);
  };
protected:
  randEngine mtRandom;
};


/**
 * \brief Abstract class for annealed criteria.
 */
class AnnealedCriteria: public Criteria
{
public:
  AnnealedCriteria(NonParametricProcess *proc):
    Criteria(proc), nCalls(0) {};
  virtual ~AnnealedCriteria(){};
  virtual void resetAnneal() { nCalls = 0; };

protected:
  unsigned int nCalls;
};


/**
 * \brief Expected improvement criterion using Schonlau annealing.
 */
class AnnealedExpectedImprovement: public AnnealedCriteria
{
public:
  AnnealedExpectedImprovement(NonParametricProcess *proc):
    AnnealedCriteria(proc), mExp(10) {};

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

protected:
  size_t mExp;
};


/**
 * \brief Lower (upper) confidence bound using Srinivas annealing
 */
class AnnealedLowerConfindenceBound: public AnnealedCriteria
{
public:
  AnnealedLowerConfindenceBound(NonParametricProcess *proc):
    AnnealedCriteria(proc), mCoef(5.0) {};

  virtual ~AnnealedLowerConfindenceBound(){};

  inline void setBetaCoef(double betac) { mCoef = betac; };

  double operator()( const vectord &x)
  {
    ++nCalls;
    size_t nDims = x.size();
    
    double beta = sqrt(2*log(static_cast<double>(nCalls*nCalls))*(nDims+1) 
		       + log(static_cast<double>(nDims))*nDims*mCoef);

    return mProc->lowerConfidenceBound(x,beta);
  };
protected:
  double mCoef;
};

////////////////////////////////////////////////////////////////////////////
/**
 * \brief Abstract class for Metacriteria interface.
 * The Metacriteria can be used to combine different criteria.
 */
class MetaCriteria
{
public:
  MetaCriteria(NonParametricProcess* proc):
    mProc(proc), mCurrentCriterium(NULL)  
  {
    criteriumSet = false;
  };

  virtual ~MetaCriteria()
  { if (criteriumSet) delete mCurrentCriterium; }

  int setCriterium(criterium_name name)
  { 
    if (criteriumSet) delete mCurrentCriterium;
    criteriumSet = false;

    mCurrentCriterium = Criteria::create(name,mProc);
    if (mCurrentCriterium == NULL) return -1;
    
    criteriumSet = true;
    return 0;  
  };
  
  inline double operator()( const vectord &x)
  { return (*mCurrentCriterium)(x); };

  virtual void reset(){};
  virtual int initializeSearch()
  {
    if ((mCurrentCriterium == NULL) || !criteriumSet)
      {
	std::cout << "Criterium not set properly" << std::endl;
	return -1;
      }
    return 0;
  };

  virtual bool checkIfBest(vectord& xNext,
			   criterium_name& name,
			   int& error_code)
  { error_code = 0;  return true; };

protected:
  NonParametricProcess* mProc;
  Criteria* mCurrentCriterium;
  bool criteriumSet;
};


/** 
 * \brief Softmax function
 * 
 * @param g gain function
 * @param eta smoothness coefficient
 * @return 
 */
inline double softmax(double g, double eta) {return exp(eta*g);};


/**
 * \brief GP_Hedge model as describen in Hoffman et al.
 *
 * The coefficients of the bandit algorithm has been carefully selected
 * according to Shapire et al. Also, the implementation has been made to
 * avoid over or underflow.
 */
class GP_Hedge: public MetaCriteria
{
public:
  GP_Hedge(NonParametricProcess *proc);
  virtual ~GP_Hedge();

  inline void reset()
  { gain = zvectord(N_ALGORITHMS_IN_GP_HEDGE); };

  int update_hedge();
  int initializeSearch()
  {
    mIndex = 0;
    mCurrentCriterium = mCriteriaList[mIndex];
    mBestLists.clear();
    return 1;
  };

  bool checkIfBest(vectord& best,criterium_name& name,int& error_code);
 
protected:
  virtual double computeLoss(const vectord& query)
  {	
    double mean, std;
    mProc->prediction(query,mean,std);
    return mean;
  }

protected:
  randEngine mtRandom;
  randFloat sampleUniform;
  vectord loss, gain, prob, cumprob;
  std::vector<Criteria*> mCriteriaList;
  std::vector<vectord> mBestLists;
  size_t mIndex;
};

/**
 * \brief Modification of the GP_Hedge algorithm where the bandit gains are
 * random outcomes.
 */
class GP_Hedge_Random: public GP_Hedge
{
public:
  GP_Hedge_Random(NonParametricProcess *proc):
    GP_Hedge(proc) {};

  virtual ~GP_Hedge_Random() {};

protected:  
  virtual double computeLoss(const vectord& query)
  { return mProc->sample_query(query,mtRandom); }
};

//@}

#endif
