/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#ifndef  _CRITERIA_FUNCTORS_HPP_
#define  _CRITERIA_FUNCTORS_HPP_

#include <numeric>
#include <algorithm>


class Criteria
{
public:
  Criteria(NonParametricProcess *proc)
  { 
    mProc = proc; 
    mYMin = 0,0;
  };

  virtual ~Criteria(){};

  inline void setMinimum(double mYMin) { mYMin = mProc->getValueAtMinimum();};
  virtual double operator()( const vectord &x) = 0;
  virtual void resetAnneal() {};  //dummy function

protected:
  NonParametricProcess *mProc;
  double mYMin;
};



////////////////////////////////////////////////////////////////////////////////
class ExpectedImprovement: public Criteria
{
public:
  ExpectedImprovement(NonParametricProcess *proc): Criteria(proc){};

  virtual ~ExpectedImprovement(){};
  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return mProc->negativeExpectedImprovement(yPred,sPred,mYMin);
  };
};


////////////////////////////////////////////////////////////////////////////////
class LowerConfidenceBound: public Criteria
{
public:
  LowerConfidenceBound(NonParametricProcess *proc): Criteria(proc) { mBeta = 1; };

  virtual ~LowerConfidenceBound(){};
  inline void setBeta(double beta) { mBeta = beta; };

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return mProc->LowerConfidenceBound(yPred,sPred,mBeta);
  };
protected:
  double mBeta;
};

////////////////////////////////////////////////////////////////////////////////
class ProbabilityOfImprovement: public Criteria
{
public:
  ProbabilityOfImprovement(NonParametricProcess *proc): Criteria(proc) {mEpsilon = 0.01;};
  virtual ~ProbabilityOfImprovement(){};
  inline void setEpsilon(double eps) {mEpsilon = eps;};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return mProc->negativeProbabilityOfImprovement(yPred,sPred,mYMin,
						mEpsilon);
  };
protected:
  double mEpsilon;
};

////////////////////////////////////////////////////////////////////////////////
class GreedyAOptimality: public Criteria
{
public:
  GreedyAOptimality(NonParametricProcess *proc): Criteria(proc){};
  virtual ~GreedyAOptimality(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return sPred;
  };
};

////////////////////////////////////////////////////////////////////////////////
class ExpectedReturn: public Criteria
{
public:
  ExpectedReturn(NonParametricProcess *proc): Criteria(proc){};

  virtual ~ExpectedReturn(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return yPred;
  };
};

////////////////////////////////////////////////////////////////////////////////
class OptimisticSampling: public Criteria
{
public:
  OptimisticSampling(NonParametricProcess *proc): mtRandom(100u), Criteria(proc){};

  virtual ~OptimisticSampling(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred, yStar = mProc->sample_query(x,mtRandom);
    mProc->prediction(query,yPred,sPred);
    return std::min(yPred,yStar);
  };
};

////////////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////////////
class AnnealedExpectedImprovement: public AnnealedCriteria
{
public:
  AnnealedExpectedImprovement(NonParametricProcess *proc):
    AnnealedCriteria(proc), mExp(10) {};

  virtual ~AnnealedExpectedImprovement(){};

  inline void setExponent(double exp) {mExp = exp;};
  void resetAnneal() { nCalls = 0; mExp(10);};

  double operator()( const vectord &x)
  {
    ++nCalls;

    if (nCalls % 10)
      mExp = std::min(1,static_cast<int>(round(mExp/2.0)));

    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return mProc->negativeExpectedImprovement(yPred,sPred,mYMin,mExp);
  };

protected:
  double mExp;
};

////////////////////////////////////////////////////////////////////////////////
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
		       + log(static_cast<double>(nDims))*nDims*coef);

    double yPred, sPred;
    mProc->prediction(query,yPred,sPred);
    return mProc->LowerConfidenceBound(yPred,sPred,beta);
  };
protected:
  double mCoef;
};

////////////////////////////////////////////////////////////////////////////////

class MetaCriteria
{
public:
  MetaCriteria(BayesOptBase *opt) 
  {
    mOpt = opt;
    mCurrentCriterium = NULL;
  };

  virtual ~MetaCriteria()
  {
    if (mCurrentCriterium != NULL)
      delete mCurrentCriterium;
  }

  int setCriterium(criterium_name name)
  {
    yieldCriteria(name,mCurrentCriterium);
  }
  
  int yieldCriteria(criterium_name name,
		    Criteria *crit)
  {
    NonParametricProcess* proc = mOpt->getSurrogateFunctionPointer();

    if (crit != NULL)
      delete crit;

    switch(name)
      {
      case C_EI:     crit = new ExpectedImprovement(proc);             break;
      case C_EI_A:   crit = new AnnealedExpectedImprovement(proc);     break;
      case C_LCB:    crit = new LowerConfidenceBound(proc);            break;
      case C_LCB_A:  crit = new AnnealedLowerConfindenceBound(proc);   break;
      case C_POI:    crit = new ProbabilityOfImprovement(proc);        break;
      case C_GREEDY_A_OPTIMALITY: crit = new GreedyAOptimality(proc);  break;
      case C_EXPECTED_RETURN:     crit = new ExpectedReturn(proc);     break;
      case C_OPTIMISTIC_SAMPLING: crit = new OptimisticSampling(proc); break;
      case C_GP_HEDGE:            crit = new GP_Hedge(proc);           break;
      case C_ERROR:
      default:
	std::cout << "Error in criterium" << std::endl; return -1;
      }
    return 1;
  };

  inline double operator()( const vectord &x)
  {
    return (*mCurrentCriterium)(x);
  }  

  virtual int returnBest(vectord& best)
  { 
    return mOpt->findOptimal(Xnext);
  }

protected:
  BayesOptBase *mOpt;
  Criteria *mCurrentCriterium;
};


class GP_Hedge: public MetaCriteria
{
public:
  GP_Hedge(BayesOptBase *opt):
    MetaCriteria(proc), mtRandom(100u),
    sampleUniform( mtRandom, realUniformDist(0,1)),
    gain(zvectord(N_ALGORITHMS_IN_GP_HEDGE)), 
    prob(zvectord(N_ALGORITHMS_IN_GP_HEDGE)),
    cumprob(zvectord(N_ALGORITHMS_IN_GP_HEDGE))
  {
    for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
      {
	Criteria *newC = NULL;
	yieldCriteria(ALGORITHMS_IN_GP_HEDGE[i],false,newC);
	mCriteriaList.push_back(newC);
      }
  };

  virtual ~GP_Hedge()
  {
    for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
      {
	delete mCriteriaList[i];
      }
  };

  inline void resetHedgeValues()
  { gain = zvectord(N_ALGORITHMS_IN_GP_HEDGE); };

  inline double softmax(double g) {return exp(mEta*g);};

  int update_hedge()
  {
    double max_g = std::max_element(gains);
    double min_g = std::min_element(gains);

    double max_l = std::max_element(loss);

    // We just care about the differences
    // std::transform(loss.begin(), loss.end(), loss.begin(),
    // 		   std::bind2nd(std::plus<double>(), max_l));       
    loss += svectord(loss.size(),max_l)

    // Optimal eta according to Shapire
    mEta = sqrt(2*log(3)/max_g);

    // To avoid overflow
    // std::transform(gain.begin(), gain.end(), gain.begin(),
    //		   std::bind2nd(std::minus<double>(), min_g));       
    gain -= svectord(gain.size(),min_g);

    std::transform(gain.begin(), gain.end(), prob.begin(),
		   softmax());       

    //Normalize
    sum_p =std::accumulate(prob.begin(),prob.end(),0);
    //    std::transform(prob.begin(), prob.end(), prob.begin(),
    //		   std::bind2nd(std::divides<double>(), sum_p));       
    prob /= sum_p;

    //Update bandits gains
    //    std::transform(gain.begin(), gain.end(), loss.begin(), gain.begin(),
    //		   std::minus<double>());       
    gain -= loss;

    std::partial_sum(prob.begin(), prob.end(), cumprob.begin(), 
		     std::plus<double>());

    double u = sampleUniform();

    for (size_t i=0; i < cumprob.size(); ++i)
      {
	if (u < cumprob(i))
	  return i;
      }
    return -1;
  };

  int returnBest(vectord& best)
  { 
    vectord bestList(N_ALGORITHMS_IN_GP_HEDGE);
    double foo;

    for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
      {
	mCurrentCriterium = mCriteriaList[i];
	mOpt->findOptimal(bestList(i));

	// The reward of the bandit problem is the estimated outcome mean 
	// at each optimal point.
	NonParametricProcess* gp = mOpt->getSurrogateFunctionPointer();
	gp->prediction(bestList(i),loss(i),foo);
      }

    int optIndex = update_hedge();
    plotResult(optIndex);

    if (optIndex < 0)
      return optIndex;

    best = bestList(optIndex);
    return 1;
  };


  int plotResult(int optIndex)
  {
    if(mParameters.verbose_level > 0)
      {
	criterium_name name;
	if (optIndex < 0)
	  name = C_ERROR;
	else
	  name = ALGORITHMS_IN_GP_HEDGE[optIndex];
	
	mOutput << crit2str(name) << " was selected." << std::endl;
      }
  };


protected:
  randEngine mtRandom;
  randFloat sampleUniform;
  vectord loss, gains, prob, cumprob;
  std::vectord<Criteria*> mCriteriaList;
};

#endif
