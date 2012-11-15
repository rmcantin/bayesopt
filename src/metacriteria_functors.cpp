#include <numeric>
#include <algorithm>
#include "boost/bind.hpp"
#include "metacriteria_functors.hpp"
#include "log.hpp"

/** 
 * \brief Softmax function
 * 
 * @param g gain function
 * @param eta smoothness coefficient
 * @return 
 */
inline double softmax(double g, double eta) 
{return exp(eta*g);};

SingleCriteria::SingleCriteria(criterium_name name, NonParametricProcess* proc): 
  MetaCriteria(proc)  
{
  mCriterium.reset(Criteria::create(name,mProc));
  if (mCriterium == NULL) 
    {
      FILE_LOG(logERROR) << "ERROR setting the criteria";
    }
};

MetaCriteria* MetaCriteria::create(criterium_name name,
				   NonParametricProcess* proc)
{
  switch(name)
    {
    case C_EI:     
    case C_EI_A:   
    case C_LCB:    
    case C_LCB_A:  
    case C_POI:    
    case C_GREEDY_A_OPTIMALITY:
    case C_EXPECTED_RETURN:
    case C_OPTIMISTIC_SAMPLING: 
      return new SingleCriteria(name,proc); 
    
    case C_GP_HEDGE: return new GP_Hedge(proc);
    case C_GP_HEDGE_RANDOM: return new GP_Hedge_Random(proc);
    default:
      FILE_LOG(logERROR) << "Error in criterium";
      return NULL;
    }
};



GP_Hedge::GP_Hedge(NonParametricProcess *proc):
  MetaCriteria(proc), mtRandom(100u),
  sampleUniform( mtRandom, realUniformDist(0,1)),
  loss_(zvectord(N_ALGORITHMS_IN_GP_HEDGE)), 
  gain_(zvectord(N_ALGORITHMS_IN_GP_HEDGE)), 
  prob_(zvectord(N_ALGORITHMS_IN_GP_HEDGE)),
  cumprob_(zvectord(N_ALGORITHMS_IN_GP_HEDGE))
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
  double max_g = *std::max_element(gain_.begin(),gain_.end());
  double min_g = *std::min_element(gain_.begin(),gain_.end());
  double max_l = *std::max_element(loss_.begin(),loss_.end());

  // We just care about the differences
  loss_ += svectord(loss_.size(),max_l);

  // To avoid overflow
  if (std::abs(max_g) > std::abs(min_g))
    gain_ -= svectord(gain_.size(),max_g);
  else
    gain_ -= svectord(gain_.size(),min_g);

  // Optimal eta according to Shapire
  max_g = *std::max_element(gain_.begin(),gain_.end());
  double eta = (std::min)(10.0,sqrt(2*log(3)/max_g));
  std::transform(gain_.begin(), gain_.end(), prob_.begin(),
		 boost::bind(softmax,_1,eta));       
    
  //Normalize
  double sum_p =std::accumulate(prob_.begin(),prob_.end(),0);
  prob_ /= sum_p;

  //Update bandits gain
  gain_ -= loss_;

  std::partial_sum(prob_.begin(), prob_.end(), cumprob_.begin(), 
		   std::plus<double>());

  double u = sampleUniform();

  for (size_t i=0; i < cumprob_.size(); ++i)
    {
      if (u < cumprob_(i))
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
      loss_(mIndex) = computeLoss(best);
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
