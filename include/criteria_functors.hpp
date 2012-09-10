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
  inline setMinimum(double mYMin) { mYMin = gp->getValueAtMinimum();};
  virtual double operator()( const vectord &x) = 0;

protected:
  NonParametricProcess *mProc;
  double mYMin;
};



////////////////////////////////////////////////////////////////////////////////
class ExpectedImprovement: public Criteria
{
public:
  ExpectedImprovement(NonParametricProcess *proc):
    Criteria(proc){};

  virtual ~ExpectedImprovement(){};
  double operator()( const vectord &x)
  {
    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return gp->negativeExpectedImprovement(yPred,sPred,mYMin);
  };
};

////////////////////////////////////////////////////////////////////////////////
class LowerConfidenceBound: public Criteria
{
public:
  LowerConfidenceBound(NonParametricProcess *proc):
    Criteria(proc) { mBeta = 1; };

  virtual ~LowerConfidenceBound(){};
  inline void setBeta(double beta) { mBeta = beta; };

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return gp->LowerConfidenceBound(yPred,sPred,mBeta);
  };
protected:
  double mBeta;
};

////////////////////////////////////////////////////////////////////////////////
class ProbabilityOfImprovement: public Criteria
{
public:
  ProbabilityOfImprovement(NonParametricProcess *proc):
    Criteria(proc) {mEpsilon = 0.01;};
  virtual ~ProbabilityOfImprovement(){};
  inline void setEpsilon(double eps) {mEpsilon = eps;};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return gp->negativeProbabilityOfImprovement(yPred,sPred,mYMin,
						mEpsilon);
  };
protected:
  double mEpsilon;
};

////////////////////////////////////////////////////////////////////////////////
class GreedyAOptimality: public Criteria
{
public:
  GreedyAOptimality(NonParametricProcess *proc):
    Criteria(proc){};
  virtual ~GreedyAOptimality(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return sPred;
  };
};

////////////////////////////////////////////////////////////////////////////////
class ExpectedReturn: public Criteria
{
public:
  ExpectedReturn(NonParametricProcess *proc):
    Criteria(proc){};
  virtual ~ExpectedReturn(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return yPred;
  };
};

////////////////////////////////////////////////////////////////////////////////
class OptimisticSampling: public Criteria
{
public:
  OptimisticSampling(NonParametricProcess *proc):
    mtRandom(100u), Criteria(proc){};
  virtual ~OptimisticSampling(){};

  double operator()( const vectord &x)
  {
    double yPred, sPred, yStar = gp->sample_query(x,mtRandom);
    gp->prediction(query,yPred,sPred);
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

  inline void resetAnneal()
  { nCalls = 0; };

protected:
  unsigned int nCalls;
};


////////////////////////////////////////////////////////////////////////////////
class AnnealedExpectedImprovement: public AnnealedCriteria
{
public:
  AnnealedExpectedImprovement(NonParametricProcess *proc):
    AnnealedCriteria(proc), mExp(1) {};

  virtual ~AnnealedExpectedImprovement(){};

  inline void setExponent(double exp) {mExp = exp;};

  double operator()( const vectord &x)
  {
    ++nCalls;

    if (nCalls % 10)
      mExp = std::min(1,static_cast<int>(round(mExp/2.0)));

    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return gp->negativeExpectedImprovement(yPred,sPred,mYMin,mExp);
  };

protected:
  double mExp;
};

////////////////////////////////////////////////////////////////////////////////
class AnnealedLowerConfindenceBound: public AnnealedCriteria
{
public:
  AnnealedLowerConfindenceBound(NonParametricProcess *proc):
    AnnealedCriteria(proc), mBeta(1) {};

  virtual ~AnnealedLowerConfindenceBound(){};

  inline void setBeta(double beta) { mBeta = beta; };

  double operator()( const vectord &x)
  {
    ++nCalls;
    size_t nDims = x.size();
    double coef = 5.0;

    beta = sqrt(2*log(nCalls*nCalls)*(nDims+1) + log(nDims)*nDims*coef);

    double yPred, sPred;
    gp->prediction(query,yPred,sPred);
    return gp->LowerConfidenceBound(yPred,sPred,mBeta);
  };
protected:
  double mBeta;
};

////////////////////////////////////////////////////////////////////////////////
class GP_Hedge: public Criteria
{
public:
  GP_Hedge(NonParametricProcess *proc):
    AnnealedCriteria(proc), mtRandom(100u),
    sampleUniform( mtRandom, realUniformDist(0,1)),
    gain(zvectord(nAlgorithmsInGPHedge)), 
    prob(zvectord(nAlgorithmsInGPHedge)),
    cumprob(zvectord(nAlgorithmsInGPHedge))
  {};

  virtual ~GP_Hedge(){};

  inline void resetHedgeValues()
  { gain = zvectord(gain.size()); }


  inline double softmax(double g) {return exp(mEta*g);};

  criterium_name update_hedge(vectord& loss)
  {
    double max_g = std::max_element(gains);
    double min_g = std::min_element(gains);

    double max_l = std::max_element(loss);

    // We just care about the differences
    std::transform(loss.begin(), loss.end(), loss.begin(),
		   std::bind2nd(std::plus<double>(), max_l));       

    // Optimal eta according to Shapire
    mEta = sqrt(2*log(3)/max_g);

    // To avoid overflow
    std::transform(gain.begin(), gain.end(), gain.begin(),
		   std::bind2nd(std::minus<double>(), min_g));       

    std::transform(gain.begin(), gain.end(), prob.begin(),
		   softmax());       

    //Normalize
    sum_p =std::accumulate(prob.begin(),prob.end(),0);
    std::transform(prob.begin(), prob.end(), prob.begin(),
		   std::bind2nd(std::divides<double>(), sum_p));       

    //Update bandits gains
    std::transform(gain.begin(), gain.end(), loss.begin(), gain.begin(),
		   std::minus<double>());       

    std::partial_sum(prob.begin(), prob.end(), cumprob.begin(), 
		     std::plus<double>());

    double u = sampleUniform();

    for (size_t i=0; i < cumprob.size(); ++i)
      {
	if (u < cumprob(i))
	  return algorithmsInGPHedge[i];
      }
    return c_error;
  }


  double operator()( const vectord &x)
  {
  };

protected:
  randEngine mtRandom;
  randFloat sampleUniform;
  vectord gains, prob, cumprob;
};

#endif
