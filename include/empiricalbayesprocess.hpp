/** \file empiricalbayesprocess.hpp
    \brief Implementes a empirical Bayesian nonparametric process with a 
    ML, MAP or similar estimate of kernel parameters. */
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


#ifndef  _EMPIRICAL_BAYES_PROCESS_HPP_
#define  _EMPIRICAL_BAYES_PROCESS_HPP_

#include "kernelregressor.hpp"
#include "inneroptimization.hpp"

namespace bayesopt
{

  /** \addtogroup  NonParametricProcesses */
  /**@{*/


  /**
   * \brief Empirical Bayesian NonParametric process.
   */
  class EmpiricalBayesProcess: public KernelRegressor, RBOptimizable
  {
  public:
    EmpiricalBayesProcess(size_t dim, bopt_params parameters,
			  Dataset& data);
    virtual ~EmpiricalBayesProcess();

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    virtual ProbabilityDistribution* prediction(const vectord &query) = 0;
		 		 
    /** 
     * \brief Updates the kernel parameters acording with a point
     * estimate (ML, MAP, etc.)
     */
    void updateKernelParameters();

    /** 
     * \brief Computes the score (eg:likelihood) of the kernel
     * parameters.  
     * Warning: To evaluate the score, it is necessary to change the parameters
     * @param x set of parameters.  
     * @return score
     */
    double evaluate(const vectord &x);

    /** 
     * \brief Computes the score (eg:likelihood) of the current kernel
     * parameters.
     * @param query set of parameters.
     * @return score
     */
    double evaluateKernelParams();


  protected:
    /** 
     * \brief Computes the negative log likelihood of the data for all
     * the parameters.
     * @return value negative log likelihood
     */
    virtual double negativeTotalLogLikelihood() = 0;


    /** 
     * \brief Computes the negative log likelihood of the data for the
     * kernel hyperparameters.
     * @return value negative log likelihood
     */
    virtual double negativeLogLikelihood() = 0;

  private:
    /**
     * Computes the negative score of the data using cross validation.
     * @return negative score
     */
    double negativeCrossValidation();

  private:
    NLOPT_Optimization* kOptimizer;
  };



  inline double EmpiricalBayesProcess::evaluate(const vectord& x)
  { 
    mKernel.setHyperParameters(x);
    return evaluateKernelParams();
  };


  /**@}*/
  
} //namespace bayesopt

#endif

