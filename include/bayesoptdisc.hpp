/** -*- c++ -*- \file bayesoptdisc.hpp \brief Discrete Bayesian optimization */
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

/**
 * \brief Sequential Kriging Optimization using different non-parametric 
 * processes as surrogate (kriging) functions. 
 */
class BAYESOPT_API BayesOptDiscrete : public BayesOptBase
{
 public:

  /** 
   * Constructor
   * @param validSet  Set of potential inputs
   */
  BayesOptDiscrete(const vecOfvec &validSet );

  /** 
   * Constructor
   * @param validSet  Set of potential inputs
   * @param params set of parameters (see parameters.h)
   */
  BayesOptDiscrete( const vecOfvec &validSet, 
		    bopt_params params);
  
  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~BayesOptDiscrete();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the [0,1] hypercube, as a 
   * normalized representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector defined in 
   * the hypercube [0,1], it might also be used as an initial point
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  int optimize( vectord &bestPoint);
  

protected:


  /** 
   * Print data for every step according to the verbose level
   * 
   * @param iteration 
   * @param xNext 
   * @param yNext 
   * 
   * @return error code
   */
  int plotStepData(size_t iteration, const vectord& xNext,
		   double yNext);

  /** Sample a set of points to initialize GP fit
   * Use pure random sampling or uniform Latin Hypercube sampling
   * as appeared in Jones 
   * @return error code
   */
  int sampleInitialPoints();

  int findOptimal(vectord &xOpt);

protected:

  vecOfvec mInputSet;               ///< List of input points

};

/**@}*/


#endif
