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
#include "criteria_functors.hpp"
#include "log.hpp"

Criteria* Criteria::create(criterium_name name,
			   NonParametricProcess* proc)
{
  switch (name)
    {
    case C_EI:     return new ExpectedImprovement(proc);
    case C_EI_A:   return new AnnealedExpectedImprovement(proc);
    case C_LCB:    return new LowerConfidenceBound(proc);
    case C_LCB_A:  return new AnnealedLowerConfindenceBound(proc);
    case C_POI:    return new ProbabilityOfImprovement(proc);
    case C_GREEDY_A_OPTIMALITY: return new GreedyAOptimality(proc);
    case C_EXPECTED_RETURN:     return new ExpectedReturn(proc);
    case C_OPTIMISTIC_SAMPLING: return new OptimisticSampling(proc);
    default:
      FILE_LOG(logERROR) << "Error in criterium";
      return NULL;
    }
};


