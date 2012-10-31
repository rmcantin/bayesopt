#include "criteria_functors.hpp"
#include "log.hpp"

Criteria* Criteria::create(criterium_name name,
			   NonParametricProcess* proc)
{
  switch(name)
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


