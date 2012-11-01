/** \file gaussian_process.hpp 
    \brief Standard zero mean gaussian process with noisy observations */
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

#ifndef  _GAUSSIAN_PROCESS_HPP_
#define  _GAUSSIAN_PROCESS_HPP_

#include "nonparametricprocess.hpp"

 
/** \addtogroup NonParametricProcesses */
/**@{*/

/**
 * \brief Standard zero mean gaussian process with noisy observations.
 */
class GaussianProcess: public NonParametricProcess
{
public:
  GaussianProcess( double noise = DEFAULT_NOISE);
  virtual ~GaussianProcess();

  /** 
   * Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query in the hypercube [0,1] to evaluate the Gaussian process
   * @param yPred mean of the predicted Gaussian distribution
   * @param sPred std of the predicted Gaussian distribution
   * 
   * @return error code.
   */	
  int prediction(const vectord &query,
  		 double& yPred, double& sPred);

  /** 
   * \brief Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query in the hypercube [0,1] to evaluate the Gaussian process
   * @return pointer to the probability distribution.
   */	
  ProbabilityDistribution* prediction(const vectord &query);


  /** 
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  double negativeLogLikelihood();

   
  /** 
   * Expected Improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * @param g exponent (used for annealing)
   * 
   * @return negative value of the expected improvement
   */
  double negativeExpectedImprovement(const vectord &query,
				     size_t g = 1);
  

  /** 
   * Lower confindence bound. Can be seen as the inverse of the Upper 
   * confidence bound
   *
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param beta std coefficient (used for annealing)
   * 
   * @return value of the lower confidence bound
   */
  double lowerConfidenceBound(const vectord &query,
			      double beta = 1);

  /** 
   * Probability of improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * @param epsilon minimum improvement margin
   * 
   * @return negative value of the probability of improvement
   */
  double negativeProbabilityOfImprovement(const vectord &query,
					  double epsilon = 0.1);

  /** 
   * Sample outcome acording to the marginal distribution at the query point.
   * 
   * @param query query point
   * @param eng boost.random engine
   * 
   * @return outcome
   */
  double sample_query(const vectord& query, randEngine& eng);

private:
  int precomputePrediction();
  vectord mAlphaV;
};

/**@}*/

#endif
