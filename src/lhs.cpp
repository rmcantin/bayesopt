#include "lhs.hpp"

void randomPerms(int arr[], size_t size, 
		 randEngine& mtRandom)
{
  /** \brief Generates an array of ramdom permutations
   * It is used to generate a uniform Latin hypercube
   */
  
  int tempIndex, randIndex;  
  randInt sample(mtRandom, intUniformDist(0,size-1));
			
  for (size_t i = 0; i < size; i++)
    arr[i] = i + 1;
  
  for (size_t last = size-1; last > 0; last--)
    {
      randIndex = sample();
      tempIndex = arr[randIndex];
      arr[randIndex] = arr[last];
      arr[last] = tempIndex;
    }
} // randomPerms 
