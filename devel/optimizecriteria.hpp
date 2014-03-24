/** \file optimizecriteria.hpp 
    \brief Class to continuous optimize criteria parameters */
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

#ifndef __OPTIMIZECRITERIA_HPP__
#define __OPTIMIZECRITERIA_HPP__

#include "inneroptimization.hpp"
#include "bayesopt.hpp"

namespace bayesopt {

  class OptimizeCriteria: public NLOPT_Optimization
  {
  public:
    explicit OptimizeCriteria(Criteria* crit):
      NLOPT_Optimization(), mCrit(crit) {};
    virtual ~OptimizeCriteria(){};

    double evaluate(const vectord& query)
    {  return (*mCrit)(query);  }
    
  private:
    Criteria* mCrit;
  };

  class OptimizeCriteriaRestricted: public NLOPT_Optimization
  {
  public:
    explicit OptimizeCriteriaRestricted(ContinuousModel* model):
      NLOPT_Optimization(), model_(model) {};
    virtual ~OptimizeCriteriaRestricted(){};

    double evaluate(const vectord& query)
    {  return model_->evaluateCriteria(query);  }
    
  private:
    ContinuousModel* model_;
  };


}


#endif
