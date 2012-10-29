#include <numeric>
#include <algorithm>
#include "boost/bind.hpp"
#include "criteria_functors.hpp"

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
    case C_GP_HEDGE:
    case C_ERROR:
    default:
      std::cout << "Error in criterium" << std::endl; 
      return NULL;
    }
};



GP_Hedge::GP_Hedge(NonParametricProcess *proc):
  MetaCriteria(proc), mtRandom(100u),
  sampleUniform( mtRandom, realUniformDist(0,1)),
  loss(zvectord(N_ALGORITHMS_IN_GP_HEDGE)), 
  gain(zvectord(N_ALGORITHMS_IN_GP_HEDGE)), 
  prob(zvectord(N_ALGORITHMS_IN_GP_HEDGE)),
  cumprob(zvectord(N_ALGORITHMS_IN_GP_HEDGE))
{
  for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
    {
      mCriteriaList.push_back(Criteria::create(ALGORITHMS_IN_GP_HEDGE[i],
					       proc));
    }
};


GP_Hedge::~GP_Hedge()
{
  for(size_t i = 0; i<N_ALGORITHMS_IN_GP_HEDGE; ++i)
    {
      delete mCriteriaList[i];
    }
};

int GP_Hedge::update_hedge()
{
  double max_g = *std::max_element(gain.begin(),gain.end());
  double min_g = *std::min_element(gain.begin(),gain.end());
  double max_l = *std::max_element(loss.begin(),loss.end());

  // We just care about the differences
  loss += svectord(loss.size(),max_l);

  // To avoid overflow
  if (std::abs(max_g) > std::abs(min_g))
    gain -= svectord(gain.size(),max_g);
  else
    gain -= svectord(gain.size(),min_g);

  // Optimal eta according to Shapire
  max_g = *std::max_element(gain.begin(),gain.end());
  double eta = std::min(10.0,sqrt(2*log(3)/max_g));
  std::transform(gain.begin(), gain.end(), prob.begin(),
		 boost::bind(softmax,_1,eta));       
    
  //Normalize
  double sum_p =std::accumulate(prob.begin(),prob.end(),0);
  prob /= sum_p;

  //Update bandits gain
  gain -= loss;

  std::partial_sum(prob.begin(), prob.end(), cumprob.begin(), 
		   std::plus<double>());

  double u = sampleUniform();

  for (size_t i=0; i < cumprob.size(); ++i)
    {
      if (u < cumprob(i))
	return i;
    }
  return -1;
};


bool GP_Hedge::checkIfBest(vectord& best, 
			   criterium_name& name,
			   int& error_code)
{ 
  if (mIndex < N_ALGORITHMS_IN_GP_HEDGE)
    {
      loss(mIndex) = computeLoss(best);
      mBestLists.push_back(best);
      error_code = 0;
      ++mIndex;
      if (mIndex < N_ALGORITHMS_IN_GP_HEDGE)
	mCurrentCriterium = mCriteriaList[mIndex];
      return false;
    }
  else
      {
	int optIndex = update_hedge();
	name = ALGORITHMS_IN_GP_HEDGE[optIndex];
	
	if (optIndex < 0)
	  error_code = optIndex;
	else
	  {
	    best = mBestLists[optIndex];
	    error_code = 0;
	  }
	return true;	
      }

  };
