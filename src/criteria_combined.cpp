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
#include <numeric>
#include <algorithm>
#include "boost/bind.hpp"
#include "log.hpp"
#include "criteria_combined.hpp"

namespace bayesopt
{

  /** 
   * \brief Softmax function
   * 
   * @param g gain function
   * @param eta smoothness coefficient
   * @return 
   */
  inline double softmax(double g, double eta) 
  {
    return exp(eta*g);
  };

  //////////////////////////////////////////////////////////////////////
  GP_Hedge::GP_Hedge():
    mtRandom(100u),
    sampleUniform( mtRandom, realUniformDist(0,1))
  {};

  int GP_Hedge::init(NonParametricProcess *proc, 
		     const std::vector<Criteria*>& list) 
  { 
    mProc = proc;
    mCriteriaList = list;
    size_t n = mCriteriaList.size();
    loss_ = zvectord(n); 
    gain_ = zvectord(n); 
    prob_ = zvectord(n);
    cumprob_ = zvectord(n);
    return 0; 
  };

  void GP_Hedge::reset()
  {
    mIndex = 0;
    mCurrentCriterium = mCriteriaList[mIndex];
    mBestLists.clear();
  };

  bool GP_Hedge::checkIfBest(vectord& best, 
			     std::string& name,
			     int& error_code)
  { 
    if (mIndex < mCriteriaList.size())
      {
	loss_(mIndex) = computeLoss(best);
	mBestLists.push_back(best);
	error_code = 0;
	++mIndex;
	if (mIndex < mCriteriaList.size())
	  mCurrentCriterium = mCriteriaList[mIndex];
	return false;
      }
    else
      {
	int optIndex = update_hedge();
	name = mCriteriaList[optIndex]->name();
      
	if (optIndex >= 0)
	  {
	    best = mBestLists[optIndex];
	    error_code = 0;
	  }
	else
	  {
	    error_code = optIndex; 
	  }
	return true;	
      }

  };


  int GP_Hedge::update_hedge()
  {
    double max_g = *std::max_element(gain_.begin(),gain_.end());
    double min_g = *std::min_element(gain_.begin(),gain_.end());
    double max_l = *std::max_element(loss_.begin(),loss_.end());

    // We just care about the differences
    loss_ += svectord(loss_.size(),max_l);

    // To avoid overflow
    if (std::abs(max_g) > std::abs(min_g))
      gain_ -= svectord(gain_.size(),max_g);
    else
      gain_ -= svectord(gain_.size(),min_g);

    // Optimal eta according to Shapire
    max_g = *std::max_element(gain_.begin(),gain_.end());
    double eta = (std::min)(10.0,sqrt(2.0*log(3.0)/max_g));
    std::transform(gain_.begin(), gain_.end(), prob_.begin(),
		   boost::bind(softmax,_1,eta));       
    
    //Normalize
    double sum_p =std::accumulate(prob_.begin(),prob_.end(),0);
    prob_ /= sum_p;

    //Update bandits gain
    gain_ -= loss_;

    std::partial_sum(prob_.begin(), prob_.end(), cumprob_.begin(), 
		     std::plus<double>());

    double u = sampleUniform();

    for (size_t i=0; i < cumprob_.size(); ++i)
      {
	if (u < cumprob_(i))
	  return i;
      }
    return -1;
  };


} //namespace bayesopt
