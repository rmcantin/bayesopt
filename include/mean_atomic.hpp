/** \file mean_atomic.hpp \brief Atomic (simple) parametric functions */
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
#include "mean_functors.hpp"

/**\addtogroup ParametricFunctions
 * @{
 */

/** \brief Abstract class for an atomic kernel */
class AtomicFunction : public ParametricFunction
{
public:
  virtual int init(size_t input_dim)
  {
    n_inputs = input_dim;
    return 0;
  };
  void setParameters(const vectord &theta) 
  {
    assert(theta.size() == n_params);
    mParameters = theta;
  };
  vectord getParameters() {return mParameters;};
  size_t nParameters() {return n_params;};

  virtual ~AtomicKernel(){};

protected:
  size_t n_params;
  size_t n_features;
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
