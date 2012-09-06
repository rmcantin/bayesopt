/** -*- c++ -*- \file elementwiseUblas.hpp \brief elementwise operations */
/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/
#ifndef  _ELEMENTWISE_UBLAS_HPP_
#define  _ELEMENTWISE_UBLAS_HPP_

// BOOST Libraries
#include <boost/numeric/ublas/vector.hpp>
#include <algorithm>

/*
template <class D>
D op_add(D in1, D in2) {return in1+in2;}
template <class D>
D op_prod(D in1, D in2) {return in1*in2;}
template <class D>
D op_div(D in1, D in2) {return in1/in2;}
*/

/** 
 * Computes the elementwise sum of two vectors.
 * 
 * @param a first vector
 * @param b second vector
 * 
 * @return sum vector.
 */
template <class v1, class v2>
v1 ublas_elementwise_add(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::plus<D>());
  return c;
}

template <class v1, class v2>
v1 ublas_elementwise_substract(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::minus<D>());
  return c;
}


template <class v1, class v2>
v1 ublas_elementwise_prod(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::multiplies<D>());
  return c;
}

template <class v1, class v2>
v1 ublas_elementwise_div(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::divides<D>());
  return c;
}

#endif
