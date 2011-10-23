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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution

#include "specialtypes.hpp"
#include "randgen.hpp"

using boost::math::factorial;
using boost::math::normal; // typedef provides default type is double.

enum criterium_name{
  expectedImprovement,
  lcb,
  poi,
  gp_hedge,
  greedyAOptimality,
  expectedReturn
};


class Criteria
{
public:

  Criteria():
    mtRandom(100u)
  {
    criterium = expectedImprovement;
    resetAnnealValues();
    resetHedgeValues();
    eta = 100;
    use_annealing = false;
  }

  virtual ~Criteria(){}

  
  inline void resetHedgeValues()
  {   g_ei = 0;  g_lcb = 0;  g_poi = 0; }
  
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

    switch (criterium)
      {
      case expectedImprovement: return negativeExpectedImprovement(yPred,sPred,yMin);
      case lcb: return lowerConfidenceBound(yPred,sPred);
      case poi: return negativeProbabilityOfImprovement(yPred,sPred,yMin);
      case greedyAOptimality: return sPred;
      case expectedReturn: return yPred;
      default: std::cout << "Error in criterium" << std::endl; return 0.0;
      }

  }

  /** 
   * Update the accumulated rewards for the GP-Hedge algorithm. See References
   * for a detailled explanation
   * 
   * @param r_ei   Reward of the Expected Improvement algorithm
   * @param r_lcb  Reward of the Lower Confidence Bound algorithm
   * @param r_poi  Reward of the Probability of Improvement algorithm
   * 
   * @return name of the selected algorithm.
   */
  criterium_name update_hedge(double r_ei, double r_lcb, double r_poi)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    double p_ei = exp(eta*g_ei);
    double p_lcb = exp(eta*g_lcb);
    double p_poi = exp(eta*g_poi);
    double sum_p = p_ei + p_lcb + p_poi;

    // Compute probabilities of choosing action
    p_ei /= sum_p; p_lcb /= sum_p; p_poi /= p_poi;
    
    // Update accumulated rewards for next time
    g_ei += r_ei; g_lcb += r_lcb; g_poi += g_poi;

    double u = sample();

    if (u < p_ei)
      return expectedImprovement;
    else if (u < p_lcb+p_ei)
      return lcb;
    else
    return poi;
  }

protected:

  /** 
   * Expected Improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * 
   * @return negative value of the expected improvement
   */
  inline double negativeExpectedImprovement(double yPred, double sPred,
					    double yMin)
  {
    double yDiff = yMin - yPred; 
    double yNorm = yDiff / sPred;
  
    if (g == 1)
	return -1.0 * ( yDiff * cdf(d,yNorm) + sPred * pdf(d,yNorm) );
    else
    {
      
      double pdfD = pdf(d,yNorm);
      double Tm2 = cdf(d,yNorm);
      double Tm1 = pdfD;
      double fg = factorial<double>(g);
      double Tact;
      double sumEI = pow(yNorm,g)*Tm2 - g*pow(yNorm,g-1)*Tm1;

      for (unsigned int ii = 2; ii < g; ii++) 
	{
	  Tact = (ii-1)*Tm2 - pdfD*pow(yNorm,ii-1);
	  sumEI += pow(-1.0,ii)*
	    (fg / ( factorial<double>(ii)* factorial<double>(g-ii) ) )*
	    pow(yNorm,g-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      return -1.0 * pow(sPred,g) * sumEI;
    }
         
  }  // negativeExpectedImprovement


  /** 
   * Lower confindence bound. Can be seen as the inverse of the Upper 
   * confidence bound
   *
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * 
   * @return value of the lower confidence bound
   */
  inline double lowerConfidenceBound(double yPred, double sPred)
  {    
    return yPred - beta*sPred;;
  }
  
  /** 
   * Probability of improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * 
   * @return negative value of the probability of improvement
   */
  inline double negativeProbabilityOfImprovement(double yPred, double sPred,
					 double yMin)
  {
    return -cdf(d,(yMin - yPred + epsilon)/sPred);
  }

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

  normal d;


};


#endif
