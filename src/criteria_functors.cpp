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
#include "log.hpp"
#include "parser.hpp"
#include "criteria_functors.hpp"
#include "criteria_atomic.hpp"
#include "criteria_combined.hpp"

namespace bayesopt
{
  CriteriaFactory::CriteriaFactory()
  {
    registry["cEI"] = & create_func<ExpectedImprovement>;
    registry["cBEI"] = & create_func<BiasedExpectedImprovement>;
    registry["cEIa"] = & create_func<AnnealedExpectedImprovement>;
    registry["cLCB"] = & create_func<LowerConfidenceBound>;
    registry["cLCBa"] = & create_func<AnnealedLowerConfindenceBound>;
    registry["cPOI"] = & create_func<ProbabilityOfImprovement>;
    registry["cAopt"] = & create_func<GreedyAOptimality>;
    registry["cExpReturn"] = & create_func<ExpectedReturn>;
    registry["cOptimisticSampling"] = & create_func<OptimisticSampling>;
    registry["cThompsonSampling"] = & create_func<ThompsonSampling>;
    registry["cDistance"] = & create_func<InputDistance>;

    registry["cSum"] = & create_func<SumCriteria>;
    registry["cProd"] = & create_func<ProdCriteria>;
    registry["cHedge"] = & create_func<GP_Hedge>;
    registry["cHedgeRandom"] = & create_func<GP_Hedge_Random>;
  }

  /// Factory method for criterion functions.
  // Criteria* CriteriaFactory::create(criterium_name name,
  // 				    KernelRegressor* proc)
  // {
  //   Criteria* c_ptr;
  //   std::vector<Criteria*> list;
  //   switch (name)
  //     {
  //     case C_EI:     c_ptr = new ExpectedImprovement(); break;
  //     case C_EI_A:   c_ptr = new AnnealedExpectedImprovement(); break;
  //     case C_LCB:    c_ptr = new LowerConfidenceBound(); break;
  //     case C_LCB_A:  c_ptr = new AnnealedLowerConfindenceBound(); break;
  //     case C_POI:    c_ptr = new ProbabilityOfImprovement(); break;
  //     case C_GREEDY_A_OPTIMALITY: c_ptr = new GreedyAOptimality(); break;
  //     case C_EXPECTED_RETURN: c_ptr = new ExpectedReturn(); break;
  //     case C_OPTIMISTIC_SAMPLING: c_ptr = new OptimisticSampling(); break;
  //     case C_GP_HEDGE: c_ptr = new GP_Hedge(); break;
  //     case C_GP_HEDGE_RANDOM: c_ptr = new GP_Hedge_Random(); break;
  //     default:
  // 	FILE_LOG(logERROR) << "Error in criterium";
  // 	return NULL;
  //     }
  //   if ((name = C_GP_HEDGE) || (name = C_GP_HEDGE_RANDOM))
  //     {
  // 	for(size_t i = 0; i < N_ALGORITHMS_IN_GP_HEDGE; ++i)
  // 	  {
  // 	    list.push_back(create(ALGORITHMS_IN_GP_HEDGE[i],proc)); 
  // 	  }
  // 	c_ptr->init(proc,list);
  //     }
  //   else
  //     {
  // 	c_ptr->init(proc);
  //     }
  //   return c_ptr;
  // };


  /** 
   * \brief Factory model for criterion functions
   * This function is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   * @param name string with the criteria structure
   * @param pointer to surrogate model
   * @return criteria pointer
   */
  Criteria* CriteriaFactory::create(std::string name,
				    NonParametricProcess* proc)
  {
    Criteria *cFunc;
    std::string os;
    std::vector<std::string> osc;
    utils::parseExpresion(name,os,osc);

    std::map<std::string,CriteriaFactory::create_func_definition>::iterator it = registry.find(os);
    if (it == registry.end()) 
      {
	FILE_LOG(logERROR) << "Error: Fatal error while parsing "
			   << "kernel function: " << os 
			   << " not found" << std::endl;
	return NULL;
      } 
    cFunc = it->second();
    if (osc.size() == 0) 
      {
	cFunc->init(proc);
      } 
    else 
      {
	std::vector<Criteria*> list;
	for(size_t i = 0; i < osc.size(); ++i)
	  {
	    list.push_back(create(osc[i],proc)); 
	  }
	cFunc->init(proc,list);
      }
    return cFunc;
  };


} //namespace bayesopt
