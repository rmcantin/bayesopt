/** -*- c++ -*- \file elementwise_ublas.hpp \brief Elementwise operations */
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

/** 
 * Computes the elementwise product of two vectors or matrices.
 * 
 * c_i = a_i * b_i
 *
 */
template <class v1, class v2>
v1 ublas_elementwise_prod(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::multiplies<D>());
  return c;
}

/** 
 * Computes the elementwise division of two vectors or matrices.
 * 
 * c_i = a_i / b_i
 *
 */
template <class v1, class v2>
v1 ublas_elementwise_div(const v1& a, const v2& b)
{
  typedef typename v1::value_type D;
  v1 c(a.size());
  std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::divides<D>());
  return c;
}

#endif
