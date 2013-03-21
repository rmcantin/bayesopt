/** \file gaussian_process.hpp 
    \brief Standard zero mean gaussian process with noisy observations */
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

#ifndef  _GAUSSIAN_PROCESS_HPP_
#define  _GAUSSIAN_PROCESS_HPP_

#include "gauss_distribution.hpp"
#include "nonparametricprocess.hpp"


namespace bayesopt
{
  
  /** \addtogroup NonParametricProcesses */
  /**@{*/

  /**
   * \brief Standard zero mean gaussian process with noisy observations.
   */
  class GaussianProcess: public NonParametricProcess
  {
  public:
    GaussianProcess(size_t dim, double noise);
    virtual ~GaussianProcess();

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    ProbabilityDistribution* prediction(const vectord &query);

  private:
    /** 
     * \brief Computes the negative log likelihood and its gradient of the data.
     * 
     * \f[ \log p(y|x,\theta,f) \propto \frac{1}{2} 
     *                             \left( y^T (K+\sigma I)^{-1} y + 
     *                             \log|K+\sigma I| \right) 
     * \f]
     *
     * @return value negative log likelihood
     */
    double negativeLogLikelihood();

    /** 
     * \brief Precompute some values of the prediction that do not depends on
     * the query
     * @return error code
     */
    int precomputePrediction();

  private:
    vectord mAlphaV;
    GaussianDistribution* d_;     ///< Pointer to distribution function
  };

  /**@}*/

} //namespace bayesopt
 

#endif
