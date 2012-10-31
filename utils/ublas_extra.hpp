/** -*- c++ -*- \file ublas_extra.hpp \brief Extra functions for Ublas library */
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

#ifndef __UBLAS_EXTRA_HPP__
#define __UBLAS_EXTRA_HPP__

#include <typeinfo>

template<class V, class D>
int append(V& vect, D element)
{
  typedef typename V::value_type VD;
  assert(typeid(VD) == typeid(D));

  // This method is super inefficient but there seems to be the uBlas style.
  const size_t size = vect.size();
  vect.resize(size+1,true);
  vect(size) = element;
  return 0;
};

#endif
