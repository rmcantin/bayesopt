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

#ifndef  _MEAN_FUNCTORS_HPP_
#define  _MEAN_FUNCTORS_HPP_

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "specialtypes.hpp"

class ParametricFunction
{
public:
  virtual void setParameters(const vectord& params)
  { mParameters = params; };
  virtual double operator() (const vectord& x) = 0;
  virtual ~ParametricFunction(){};

protected:
  vectord mParameters;
};

///////////////////////////////////////////////////////////////////////////

class ZeroFunction: public ParametricFunction
{
public:
  double operator() (const vectord& x)
  { return 0.0; };
};

///////////////////////////////////////////////////////////////////////////

class OneFunction: public ParametricFunction
{
public:
  double operator() (const vectord& x)
  { return 1.0; };
};


///////////////////////////////////////////////////////////////////////////

class ConstantFunction: public ParametricFunction
{
public:
  double operator() (const vectord& x)
  { return mParameters(0); };
};

///////////////////////////////////////////////////////////////////////////

class LinearFunction: public ParametricFunction
{
public:
  double operator() (const vectord& x)
  { return boost::numeric::ublas::inner_prod(x,mParameters);  };
};

///////////////////////////////////////////////////////////////////////////

class LinearPlusConstantFunction: public ParametricFunction
{
public:
  void setParameters(const vectord& params)
  { 
    mConstParam = params(0);
    mParameters = boost::numeric::ublas::project(params, 
			boost::numeric::ublas::range(1, params.size())); 
  };
  
  double operator() (const vectord& x)
  { return boost::numeric::ublas::inner_prod(x,mParameters) + mConstParam;  };

protected:
  double mConstParam;
};


#endif
