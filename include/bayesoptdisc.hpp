/**  \file bayesoptdisc.hpp \brief Discrete Bayesian optimization */
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

#ifndef  _BAYESOPTDISC_HPP_
#define  _BAYESOPTDISC_HPP_

#include "bayesoptbase.hpp"

/** \addtogroup BayesOpt */
/*@{*/

namespace bayesopt
{

  /**
   * \brief Sequential Kriging Optimization using different non-parametric 
   * processes as surrogate (kriging) functions. 
   */
  class BAYESOPT_API DiscreteModel : public BayesOptBase
  {
  public:

    /** 
     * Constructor
     * @param validSet  Set of potential inputs
     */
    DiscreteModel(const vecOfvec &validSet );

    /** 
     * Constructor
     * @param validSet  Set of potential inputs
     * @param params set of parameters (see parameters.h)
     */
    DiscreteModel( const vecOfvec &validSet, 
		 bopt_params params);
    
    /** Default destructor  */
    virtual ~DiscreteModel();

    /** Initialize the optimization process. */
    void initializeOptimization();

    /** 
     * Once the optimization has been perfomed, return the optimal
     * point.
     */
    vectord getFinalResult();

    
  protected:
    
    
    /** 
     * Print data for every step according to the verbose level
     * 
     * @param iteration 
     * @param xNext 
     * @param yNext 
     */
    void plotStepData(size_t iteration, const vectord& xNext,
		     double yNext);

    /** 
     * Sample a set of points to initialize GP fit.
     * Use pure random sampling or uniform Latin Hypercube sampling
     * as appeared in Jones 
     * @return error code
     */
    int sampleInitialPoints();

    /** 
     * \brief Wrapper for the target function normalize in the hypercube
     * [0,1]
     * @param query point to evaluate in [0,1] hypercube
     * @return actual return value of the target function
     */
    inline double evaluateSampleInternal( const vectord &query )
    { return evaluateSample(query); }; 

    int findOptimal(vectord &xOpt);

  protected:
    vecOfvec mInputSet;               ///< List of input points

  };

} //namespace bayesopt

/**@}*/


#endif
