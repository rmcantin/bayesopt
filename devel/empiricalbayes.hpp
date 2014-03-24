/** \file empiricalbayes.hpp
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

#ifndef _EMPIRICAL_BAYES_HPP_
#define _EMPIRICAL_BAYES_HPP_

namespace bayesopt
{

  /** \addtogroup  LearningMethods */
  /**@{*/


  /**
   * \brief Empirical Bayesian NonParametric process.
   */
  class ConditionalBayesProcess: public KernelRegressor, RBOptimizable
  {
  public:
  

} // namespace bayesopt

#endif
