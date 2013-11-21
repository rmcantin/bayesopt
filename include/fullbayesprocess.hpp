/** \file fullbayesprocess.hpp
    \brief Implementes a fully Bayesian nonparametric process with a 
    sampling distribution over kernel parameters. */
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


#ifndef  _FULL_BAYES_PROCESS_HPP_
#define  _FULL_BAYES_PROCESS_HPP_

namespace bayesopt
{

  /** \addtogroup  NonParametricProcesses */
  /**@{*/


  /**
   * \brief Full Bayesian NonParametric process.
   */
  class FullBayesProcess: public NonParametricProcess
  {
  public:
    static const size_t N_PROC = 10;

    FullBayesProcess(size_t dim, bopt_params params);
    virtual ~FullBayesProcess();

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    ProbabilityDistribution* prediction(const vectord &query);

    /** 
     * \brief Updates the kernel parameters acording with a point
     * estimate (ML, MAP, etc.)
     * @return error code
     */
    int updateKernelParameters();

  private:
    std::vector<NonParametricProcess*>   mVProc;
    vectord                            mWeights;
    
    MixtureDistribution* d_;      //!< Predictive distributions
  };


  /**@}*/

} //namespace bayesopt


#endif
