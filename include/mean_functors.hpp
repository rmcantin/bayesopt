/** \file mean_functors.hpp \brief Mean (parametric) functions. */
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

#ifndef  _MEAN_FUNCTORS_HPP_
#define  _MEAN_FUNCTORS_HPP_

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "parameters.h"
#include "specialtypes.hpp"

/**\addtogroup ParametricFunctions
 *  \brief Set of parametric models for surrogate modelling
 */
//@{

/** \brief Interface for mean functors */
class ParametricFunction
{
public:
  static ParametricFunction* create(mean_name name, const vectord& params);
  virtual void setParameters(const vectord& params){ mParameters = params; };
  virtual double getMean(const vectord& x) = 0;
  virtual vectord getFeatures(const vectord& x) = 0;  
  virtual size_t getN(const vectord& x) = 0;  
  
  virtual vectord operator()(const vecOfvec& x)
  {
    vectord result(x.size());
    std::vector<vectord>::const_iterator x_it = x.begin();
    std::vector<vectord>::const_iterator x_end = x.end();
    vectord::iterator res_it = result.begin();
    while(x_it != x_end)
      {
	*res_it++ = getMean(*x_it++);
      }
    return result;
  };

  virtual matrixd getAllFeatures(const vecOfvec& x)
  {
    size_t nf = getN(x[0]);
    matrixd result(nf,x.size());

    for(size_t ii=0; ii< x.size(); ++ii)
      {
	column(result,ii) = getFeatures(x[ii]);
      }

    return result;
  };
  
  virtual ~ParametricFunction(){};

protected:
  vectord mParameters;
};


/** \brief Constant zero function */
class ZeroFunction: public ParametricFunction
{
public:
  double getMean (const vectord& x) { return 0.0; };
  vectord getFeatures(const vectord& x) { return zvectord(1); };  
  size_t getN(const vectord& x) {return 1;}  
};


/** \brief Constant one function */
class OneFunction: public ParametricFunction
{
public:
  double getMean (const vectord& x) { return 1.0; };
  vectord getFeatures(const vectord& x) { return svectord(1,1.0); };  
  size_t getN(const vectord& x) {return 1;}  
};


/** \brief Constant function. 
    The first parameter indicates the constant value. */
class ConstantFunction: public ParametricFunction
{
public:
  double getMean (const vectord& x) { return mParameters(0); };
  vectord getFeatures(const vectord& x) { return svectord(1,1.0); };  
  size_t getN(const vectord& x) {return 1;}  
};


/** \brief Linear combination function. 
    Each parameter indicates the coefficient of each dimension. */
class LinearFunction: public ParametricFunction
{
public:
  double getMean (const vectord& x)
  { return boost::numeric::ublas::inner_prod(x,mParameters);  };
  vectord getFeatures(const vectord& x) { return x; };  
  size_t getN(const vectord& x) {return x.size();}  
};


/** \brief Linear combination plus constant function. 
    The first parameter indicates the constant value. */
class LinearPlusConstantFunction: public ParametricFunction
{
public:
  void setParameters(const vectord& params)
  { 
    mConstParam = params(0);
    mParameters = boost::numeric::ublas::project(params, 
			boost::numeric::ublas::range(1, params.size())); 
  };
  
  double getMean (const vectord& x)
  { return boost::numeric::ublas::inner_prod(x,mParameters) + mConstParam;  };

  vectord getFeatures(const vectord& x) 
  {
    using boost::numeric::ublas::range;
    using boost::numeric::ublas::project;
    vectord res(x.size()+1);
    res(0) = 1;
    project(res,range(1,res.size())) = x;
    return res; 
  };  

  size_t getN(const vectord& x) {return 1+x.size();}  

protected:
  double mConstParam;
};

//@}

#endif
