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

#ifndef  _CRITERIA_HPP_
#define  _CRITERIA_HPP_

#include <algorithm>

#include "ctypes.h"
#include "specialtypes.hpp"
#include "randgen.hpp"


class Criteria
{
public:

  Criteria():
    mtRandom(100u)
  {
    criterium = c_ei;
    resetAnnealValues();
    resetHedgeValues();
    eta = 1;
    use_annealing = false;
  }

  virtual ~Criteria(){}

  
  inline void resetHedgeValues()
  { g_ei = 0;  g_lcb = 0;  g_poi = 0; }
  
  inline void resetAnnealValues()
  { n_calls = 0; g = 1; beta = 1; epsilon = 0.01; }

  inline void setCriterium(criterium_name c)
  { criterium = c; }

  inline void setAnnealing(bool anneal)
  { use_annealing = anneal; }

  inline double evaluate(NonParametricProcess &gp, const vectord &query)
  {
    n_calls++;
    if (use_annealing) updateCoolingScheme(query.size());

    double yPred, sPred, yMin = gp.getValueAtMinimum(); 
    gp.prediction(query,yPred,sPred);
    double yStar;
    
     switch (criterium)
      {
      case c_ei: return gp.negativeExpectedImprovement(yPred,sPred,yMin,g);
      case c_lcb: return gp.lowerConfidenceBound(yPred,sPred,beta);
      case c_poi: return gp.negativeProbabilityOfImprovement(yPred,sPred,yMin,epsilon);
      case c_greedyAOptimality: return sPred;
      case c_expectedReturn: return yPred;
      case c_optimisticSampling: 
	yStar = gp.sample_query(query,mtRandom);
	return std::min(yPred,yStar);
      case c_gp_hedge:
      default: std::cout << "Error in criterium" << std::endl; return 0.0;
      }

  }

  inline double abs_max(double a, double b)
  { return (std::abs(a)<std::abs(b)) ? b:a; }

  /** 
   * Update the accumulated rewards for the GP-Hedge algorithm. See References
   * for a detailled explanation. Although it has been carefully implemented to 
   * avoid overflow or underflow it may still suffer from it if the predictions 
   * of the GP produce very large numbers. 
   * 
   * @param l_ei   Loss of the Expected Improvement algorithm
   * @param l_lcb  Loss of the Lower Confidence Bound algorithm
   * @param l_poi  Loss of the Probability of Improvement algorithm
   * 
   * @return name of the selected algorithm.
   */
  criterium_name update_hedge(double l_ei, double l_lcb, double l_poi)
  {
    // TODO: Vectorize the function to support an arbitrary number of criteria.
    randFloat sample( mtRandom, realUniformDist(0,1) );

    double max_g = std::max(g_ei,std::max(g_lcb,g_poi));
    double min_g = std::min(g_ei,std::min(g_lcb,g_poi));

    double max_l = std::max(l_ei,std::max(l_lcb,l_poi));
    l_ei += max_l; l_lcb += max_l; l_poi += max_l;        // We just care about the differences

    eta = sqrt(2*log(3)/max_g);                           // Optimal eta according to Shapire

    // To avoid overflow
    double offset = min_g;
    g_ei -= offset; g_lcb -= offset; g_poi -= offset; 

    double p_ei = exp(eta*g_ei);
    double p_lcb = exp(eta*g_lcb);
    double p_poi = exp(eta*g_poi);
    double sum_p = p_ei + p_lcb + p_poi;

    // Compute probabilities of choosing action
    p_ei  /= sum_p; p_lcb /= sum_p; p_poi /= sum_p;
    
    // Update accumulated rewards for next time
    g_ei -= l_ei; g_lcb -= l_lcb; g_poi -= l_poi; 

    double u = sample();

    if (u < p_ei)
      return c_ei;
    else if (u < p_lcb+p_ei)
      return c_lcb;
    else
    return c_poi;
  }

protected:

  /** 
   * Updates the parameters that can be used for annealing
   * 
   * @param ndims # of imput dimensions.
   */
  inline void updateCoolingScheme(size_t ndims)
  {
    double coef = 5;

    if (n_calls%10)
      g = std::max(1,static_cast<int>(round(g/2.0)));
    beta = sqrt(2*log(n_calls*n_calls)*(ndims+1) + log(ndims)*ndims*coef);
  }

  
protected:

  criterium_name criterium;
  randEngine mtRandom;
  unsigned int n_calls, g;

  double beta, epsilon, eta; 
  double g_poi,g_lcb,g_ei;
  bool use_annealing;

};


#endif
