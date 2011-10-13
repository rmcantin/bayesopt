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
  for (std::vector<int>::iterator it=arr.end(); it!=arr.begin(); --it)
    iter_swap(arr.begin()+sample(),it);
  
  for (std::vector<int>::iterator it=arr.begin(); it!=arr.end(); ++it)
    std::cout << " " << *it;
  
  /*
  for (size_t i = 0; i < size; i++)
    arr[i] = i + 1;
 
  for (size_t last = size-1; last > 0; last--)
    {
      randIndex = sample();
      tempIndex = arr[randIndex];
      arr[randIndex] = arr[last];
      arr[last] = tempIndex;
    }
  */
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
	  double perVal = static_cast<double>(perms[j]);
	  std::cout << perVal << " ";
	  Result(j,i) = (perVal  - sample() ) / ndA;
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
