/** \file mixture_distribution.hpp 
    \brief Mixture of gaussians probability distribution */
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


#ifndef __MIXTURE_DISTRIBUTION_HPP__
#define __MIXTURE_DISTRIBUTION_HPP__

#include "prob_distribution.hpp" 

class MixtureDistribution: public ProbabilityDistribution
{
public:
  MixtureDistribution(size_t n);
  virtual ~MixtureDistribution();

  /** 
   * \brief Probability density function
   * @param x query point
   * @return probability
   */
  double pdf(double x);

  /** 
   * \brief Expected Improvement algorithm for minimization
   * @param min  minimum value found
   * @param g exponent (used for annealing)
   *
   * @return negative value of the expected improvement
   */
  double negativeExpectedImprovement(double min, size_t g);

  /** 
   * \brief Lower confindence bound. Can be seen as the inverse of the Upper 
   * confidence bound
   * @param beta std coefficient (used for annealing)
   * @return value of the lower confidence bound
   */
  double lowerConfidenceBound(double beta);

  /** 
   * Probability of improvement algorithm for minimization
   * @param min  minimum value found
   * @param epsilon minimum improvement margin
   * 
   * @return negative value of the probability of improvement
   */
  double negativeProbabilityOfImprovement(double min,
					  double epsilon);

  /** 
   * Sample outcome acording to the marginal distribution at the query point.
   * @param eng boost.random engine
   * 
   * @return outcome
   */
  double sample_query(randEngine& eng);

  double getMean();
  double getStd();
  double getGaussianStd();

private:
  std::vector<ProbabilityDistribution*> mPD; 
  vectord mW; 
};



#endif
