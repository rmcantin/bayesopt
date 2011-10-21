#ifndef _LHS_HPP_
#define _LHS_HPP_

#include "randgen.hpp"


struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return ++current;}
} UniqueNumber;


/** \brief Generates an array of ramdom permutations
 * It is used to generate a uniform Latin hypercube
 */
void randomPerms(std::vector<int>& arr, 
		 randEngine& mtRandom)
{
  
  randInt sample(mtRandom, intUniformDist(0,arr.size()-1));
			
  generate (arr.begin(), arr.end(), UniqueNumber);
  for (std::vector<int>::iterator it=arr.begin(); it!=arr.end(); ++it)
    iter_swap(arr.begin()+sample(),it);
} // randomPerms 

/** \brief Latin hypercube sampling
 * It is used to generate a uniform Latin hypercube
 */
template<class M>
int lhs(M& Result,
	randEngine& mtRandom)
{
  randFloat sample( mtRandom, realUniformDist(0,1) );
  size_t nA = Result.size1();
  size_t nB = Result.size2();
  double ndA = static_cast<double>(nA);
  std::vector<int> perms(nA);
      
  for (size_t i = 0; i < nB; i++)
    {
      randomPerms(perms, mtRandom);
      
      for (size_t j = 0; j < nA; j++)
	{		
	  Result(j,i) = ( static_cast<double>(perms[j]) - sample() ) / ndA;
	}
    }

  return 1;
}

/** \brief Uniform hypercube sampling
 * It is used to generate a set of uniformly distributed
 * samples in hypercube
 */
template<class M>
int uniformSampling(M& Result,
		    randEngine& mtRandom)
{
  randFloat sample( mtRandom, realUniformDist(0,1) );
  size_t nA = Result.size1();
  size_t nB = Result.size2();

  // TODO: Improve with iterators
  for (size_t i = 0; i < nA; i++)
    for (size_t j = 0; j < nB; j++)
	  Result(i,j) = sample();

  return 1;
}


#endif
