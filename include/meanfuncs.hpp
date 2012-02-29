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

#ifndef  _MEANFUNCS_HPP_
#define  _MEANFUNCS_HPP_

#include <boost/numeric/ublas/vector.hpp>

namespace means
{
  using namespace boost::numeric::ublas;	

  /** 
   * Constant unit function
   * 
   * @return 1
   */
  inline double One( const vector<double> &x )
  { return 1; } 
  
  /** 
   * Constant zero function
   * 
   * @return 0
   */
  inline double Zero( const vector<double> &x)
  { return 0; } 

  /** 
   * Linear function
   * 
   * @param x variable
   * @param a coefficient
   * 
   * @return a \dot x
   */
  inline double Linear (const vector<double> &x,
			const vector<double> &a)
  { return inner_prod(x,a); }
}


#endif
