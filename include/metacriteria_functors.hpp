/**  \file metacriteria_functors.hpp 
    \brief Functions to combine different criteria */
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

#ifndef  _METACRITERIA_FUNCTORS_HPP_
#define  _METACRITERIA_FUNCTORS_HPP_

#include <boost/scoped_ptr.hpp>
#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup MetaCriteriaFunctions
   * \brief Functions that combine different criteria in a single step.
   */
  //@{
  
  /**
   * \brief Abstract class for Metacriteria interface.
   * The Metacriteria can be used to combine different criteria.
   */
  class MetaCriteria
  {
  public:
    MetaCriteria(NonParametricProcess* proc): 
      mProc(proc) {};

    virtual ~MetaCriteria() {};

    virtual bool requireComparison() = 0;
    virtual double operator()(const vectord &x) = 0;

    //Dummy functions.
    virtual int initializeSearch() { assert(false); };
    virtual bool checkIfBest(vectord& xNext,
			     std::string& name,
			     int& error_code) { assert(false); };
  
    //Factory
    static MetaCriteria* create(criterium_name name,
				NonParametricProcess* proc);

  protected:
    NonParametricProcess* mProc;
  };

  /// Abstract class for atomic criteria (only one function)
  class AtomicCriteria: public MetaCriteria
  {
  public:
    AtomicCriteria(NonParametricProcess* proc, criterium_name name);
    AtomicCriteria(NonParametricProcess* proc, std::string name):
      MetaCriteria(proc) { /* TODO: */ };
    virtual ~AtomicCriteria(){};
  protected:
    boost::scoped_ptr<Criteria> mCriterium;
  };

  /// Abstract class for combined criteria (multiple functions)
  class CombinedCriteria: public MetaCriteria
  {
  public:
    CombinedCriteria(NonParametricProcess* proc);
    CombinedCriteria(NonParametricProcess* proc, 
		     const std::vector<std::string>& name):
      MetaCriteria(proc) { /* TODO: */ };
    virtual ~CombinedCriteria();
  protected:
    std::vector<Criteria*> mCriteriaList;
  };


  /**
   * \brief Wrapper class for single criterion functions.
   */
  class SingleCriteria: public AtomicCriteria
  {
  public:
    SingleCriteria(criterium_name name, NonParametricProcess* proc):
      AtomicCriteria(proc,name) {};
    virtual ~SingleCriteria() {};
  
    bool requireComparison(){ return false; };
    double operator()(const vectord &x)  { return (*mCriterium)(x); };
  };


  /**
   * \brief Wrapper class for linear combination of criterion functions.
   */
  class SumCriteria: public CombinedCriteria
  {
  public:
    SumCriteria(NonParametricProcess* proc):
      CombinedCriteria(proc){};
    virtual ~SumCriteria() {};
  
    bool requireComparison(){ return false; };
    double operator()(const vectord &x)  
    {
      double sum = 0.0;
      for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
	{
	  sum += (*mCriteriaList[i])(x); 
	}
      return sum;
    };
  };


  /**
   * \brief GP_Hedge model as describen in Hoffman et al. \cite Hoffman2011
   *
   * The coefficients of the bandit algorithm has been carefully selected
   * according to Shapire et al. Also, the implementation has been made to
   * avoid over or underflow.
   */
  class GP_Hedge: public CombinedCriteria
  {
  public:
    GP_Hedge(NonParametricProcess *proc);
    virtual ~GP_Hedge() {};

    bool requireComparison(){ return true; };
    double operator()(const vectord &x) { return (*mCurrentCriterium)(x); };

    int initializeSearch();
    bool checkIfBest(vectord& best, std::string& name,int& error_code);

  protected:
    int update_hedge();

    randEngine mtRandom;
    randFloat sampleUniform;
    vectord loss_, gain_, prob_, cumprob_;
    Criteria* mCurrentCriterium;
    std::vector<vectord> mBestLists;
    size_t mIndex;

 
  private:
    double computeLoss(const vectord& query)
    { return mProc->prediction(query)->getMean(); }
  };

  /**
   * \brief Modification of the GP_Hedge algorithm where the bandit gains are
   * random outcomes (Thompson sampling).
   */
  class GP_Hedge_Random: public GP_Hedge
  {
  public:
    GP_Hedge_Random(NonParametricProcess *proc):
      GP_Hedge(proc) {};

    virtual ~GP_Hedge_Random() {};

  private:
    double computeLoss(const vectord& query)
    { return mProc->prediction(query)->sample_query(mtRandom); }
  };

  //@}

} //namespace bayesopt


#endif
