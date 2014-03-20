/**  \file indexvector.hpp \brief Generators for index vectors. */
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

#ifndef __INDEX_VECTOR_HPP__
#define __INDEX_VECTOR_HPP__

#include <algorithm>

namespace bayesopt
{
  namespace utils
  {
    /**
     * \brief Simple class to generate sequences of unique numbers
     */
    class CUnique 
    {
    private:
      int current;
    public:
      CUnique() {current=0;}
      int operator()() {return ++current;}
    };

    /** 
     * Generates a vector of indexes (0..n)
     * @param n vector size
     * @return index vector
     */
    inline std::vector<int> returnIndexVector(size_t n)
    {
      CUnique UniqueNumber;
      std::vector<int> arr(n);
      generate (arr.begin(), arr.end(), UniqueNumber);
      return arr;
    };

    /** 
     * Modify a vector of indexes (0..n)
     * @param arr vector
     */
    inline void modifyIndexVector(std::vector<int>& arr)
    {
      CUnique UniqueNumber;
      generate (arr.begin(), arr.end(), UniqueNumber);
    };

  } //namespace utils
} //namespace bayesopt

#endif
