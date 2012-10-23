/** -*- c++ -*- \file indexvector.hpp \brief Generators for index vectors. */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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

#ifndef _LHS_HPP_
#define _LHS_HPP_

#include "randgen.hpp"
#include "indexvector.hpp"

/** \brief Modify an array using ramdom permutations.
 *
 * It is used to generate a uniform Latin hypercube.
 * Equivalent to std::random_shuffle but using boost::random
 */
template<class D>
void randomPerms(D& arr, 
		 randEngine& mtRandom)
{
  typedef typename D::iterator iter;

  randInt sample(mtRandom, intUniformDist(0,arr.size()-1));
  for (iter it=arr.begin(); it!=arr.end(); ++it)
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
  //  std::vector<int> perms(nA);
  
  for (size_t i = 0; i < nB; i++)
    {
      std::vector<int> perms = returnIndexVector(nA);
      randomPerms(perms, mtRandom);
      //std::random_shuffle(perms.begin(),perms.end());
      
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
