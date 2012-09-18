/** -*- c++ -*- \file boundingbox.hpp \brief Functions for bounding box limits */
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
#ifndef  _BOUNDING_BOX_HPP_
#define  _BOUNDING_BOX_HPP_

// BOOST Libraries
#include <boost/numeric/ublas/vector.hpp>
#include "elementwise_ublas.hpp"

template <class V>
class BoundingBox
{
public:
  BoundingBox(const V &lbound, const V &ubound):
    mLowerBound(lbound), mRangeBound(ubound - lbound)
  {};

  virtual ~BoundingBox(){};

  inline V unnormalizeVector( const V &vin )
  {
    return ublas_elementwise_prod(vin,mRangeBound) + mLowerBound;
  };  // unnormalizeVector

  inline V normalizeVector( const V &vin )
  {
    return ublas_elementwise_div(vin - mLowerBound, mRangeBound);
  }  // normalizeVector
  
protected:
  V mLowerBound, mRangeBound; ///< Lower bound and range of the input space
};

#endif
