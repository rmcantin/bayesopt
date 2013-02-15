
/** \file kernel_functors.hpp \brief Kernel (covariance) functions */
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

#ifndef  _KERNEL_FUNCTORS_HPP_
#define  _KERNEL_FUNCTORS_HPP_

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "parameters.h"
#include "specialtypes.hpp"
#include "elementwise_ublas.hpp"

/**\addtogroup KernelFunctions
 * \brief Set of kernel or covariance functions for the nonparametric processes.
 */
//@{

/** \brief Interface for kernel functors */
class Kernel
{
public:
  static Kernel* create(std::string name, size_t imput_dim);
  static Kernel* create(kernel_name name, size_t imput_dim);

  virtual int init(size_t input_dim) {return 0;};
  virtual int init(size_t input_dim, Kernel* left, Kernel* right) {return 0;};

  virtual void setHyperParameters(const vectord &theta) = 0;
  virtual vectord getHyperParameters() = 0;
  virtual size_t nHyperParameters() = 0;

  virtual double operator()( const vectord &x1, const vectord &x2 ) = 0;
  virtual double gradient( const vectord &x1, const vectord &x2,
			   size_t component ) = 0;
  virtual ~Kernel(){};

protected:
  size_t n_inputs;
};

/////////////////////////////////////////////////////////////////////////////
/** \brief Abstract class for an atomic kernel */
class AtomicKernel : public Kernel
{
public:
  virtual int init(size_t input_dim)
  {
    n_inputs = input_dim;
    return 0;
  };
  void setHyperParameters(const vectord &theta) 
  {
    assert(params.size() == n_params);
    params = theta;
  };
  vectord getHyperParameters() {return params;};
  size_t nHyperParameters() {return n_params;};

  virtual ~AtomicKernel(){};

protected:
  size_t n_params;
  vectord params;
};

/////////////////////////////////////////////////////////////////////////////
/** \brief Abstract class for combined kernel */
class CombinedKernel : public Kernel
{
public:
  virtual int init(size_t input_dim, Kernel* left, Kernel* right)
  {
    n_inputs = input_dim;
    this->left = left;
    this->right = right;
    return 0;
  };
  void setHyperParameters(const vectord &theta) 
  {
    using boost::numeric::ublas::subrange;

    size_t n_lhs = left->nHyperParameters();
    size_t n_rhs = right->nHyperParameters();
    assert(params.size() == n_lhs + n_rhs);
    left->setHyperParameters(subrange(theta,0,n_lhs));
    right->setHyperParameters(subrange(theta,n_lhs,n_lhs+n_rhs));
  };

  vectord getHyperParameters() 
  {
    using boost::numeric::ublas::subrange;

    size_t n_lhs = left->nHyperParameters();
    size_t n_rhs = right->nHyperParameters();
    vectord par(n_lhs + n_rhs);
    subrange(par,0,n_lhs) = left->getHyperParameters();
    subrange(par,n_lhs,n_lhs+n_rhs) = right->getHyperParameters();
    return par;
  };

  size_t nHyperParameters() 
  {
    size_t n_lhs = left->nHyperParameters();
    size_t n_rhs = right->nHyperParameters();
    return n_lhs + n_rhs;
  };

  virtual ~CombinedKernel()
  {
    delete left;
    delete right;
  };

protected:
  Kernel* left;
  Kernel* right;
};

////////////////////////////////////////////////////////////////////////////////
/** \brief Sum of two kernels */
class KernelSum: public CombinedKernel
{
public:
  double operator()(const vectord &x1, const vectord &x2)
  {
    return (*left)(x1,x2) + (*right)(x1,x2);
  };

  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    return left->gradient(x1,x2,component) + right->gradient(x1,x2,component);
  };
};

////////////////////////////////////////////////////////////////////////////////
/** \brief Product of two kernels */
class KernelProd: public CombinedKernel
{
public:
  double operator()(const vectord &x1, const vectord &x2)
  {
    return (*left)(x1,x2) * (*right)(x1,x2);
  };

  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    return 0.0; //TODO: Not implemented
  };
};

////////////////////////////////////////////////////////////////////////////////
/** \brief Constant kernel.  
 * Combined with sum and product kernel can be used to scale a kernel
 * or add noise.
 */

class ConstKernel: public AtomicKernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 1;
    n_inputs = input_dim;
    return 0;
  };

  double operator()(const vectord &x1, const vectord &x2)
  {
    return params(0);
  };

  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    return 0.0;
  };
};

/////////////////////////////////////////////////////////////////////////////
/** \brief Linear kernel. */
class LinKernel: public AtomicKernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 0;
    n_inputs = input_dim;
    return 0;
  };

  double operator()(const vectord &x1, const vectord &x2)
  {
    assert(x1.size() == x2.size());
    return boost::numeric::ublas::inner_prod(x1,x2); 
  };

  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    assert(false);
  };
};



/////////////////////////////////////////////////////////////////////////////
/** \brief Abstract class for isotropic kernel functors */
class ISOkernel : public AtomicKernel
{
public:
  virtual ~ISOkernel(){};

protected:
  inline double computeScaledNorm2(const vectord &x1, const vectord &x2)
  {  
    assert(n_inputs == x1.size());
    assert(x1.size() == x2.size());
    return norm_2(x1-x2)/params(0); 
  };
};


/////////////////////////////////////////////////////////////////////////////
/** \brief Abstract class for anisotropic kernel functors. 
 * Typically ARD (Automatic Relevance Determination)
 */
class ARDkernel : public AtomicKernel
{
public:
  virtual ~ARDkernel(){};

protected:
  inline vectord computeScaledDiff(const vectord &x1, const vectord &x2)
  {
    assert(n_inputs == x1.size());
    assert(x1.size() == x2.size());
    assert(x1.size() == params.size());

    vectord xd = x1-x2;
    return ublas_elementwise_div(xd, params); 
  };
};



////////////////////////////////////////////////////////////////////////////////
/** \brief Matern kernel of 1st order */
class MaternIso1: public ISOkernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 1;
    n_inputs = input_dim;
    return 0;
  };

  double operator()(const vectord &x1, const vectord &x2)
  {
    double r = computeScaledNorm2(x1,x2);
    return exp(-r);
  };

  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    double r = computeScaledNorm2(x1,x2);
    return r*exp(-r);
  };
};


////////////////////////////////////////////////////////////////////////////////
/** \brief Matern kernel of 3rd order */
class MaternIso3: public ISOkernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 1;
    n_inputs = input_dim;
    return 0;
  };

  double operator()( const vectord &x1, const vectord &x2)
  {
    double r = sqrt(3.0) * computeScaledNorm2(x1,x2);
    double er = exp(-r);
    return (1+r)*er;
  };

  double gradient( const vectord &x1, const vectord &x2,
		   size_t component)
  {
    double r = sqrt(3.0) * computeScaledNorm2(x1,x2);
    double er = exp(-r);
    return r*r*er; 
  };
};


////////////////////////////////////////////////////////////////////////////////
/** \brief Matern kernel of 5th order */
class MaternIso5: public ISOkernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 1;
    n_inputs = input_dim;
    return 0;
  };

  double operator()( const vectord &x1, const vectord &x2)
  {
    double r = sqrt(5.0) * computeScaledNorm2(x1,x2);
    double er = exp(-r);
    return (1+r*(1+r/3))*er;
  };
  double gradient( const vectord &x1, const vectord &x2,
		   size_t component)
  {    
    double r = sqrt(5.0) * computeScaledNorm2(x1,x2);
    double er = exp(-r);
    return r*(1+r)/3*r*er; 
  };
};



////////////////////////////////////////////////////////////////////////////////
/** \brief Square exponential (Gaussian) kernel. Isotropic version. */
class SEIso: public ISOkernel
{
public:
  int init(size_t input_dim)
  {
    n_params = 1;
    n_inputs = input_dim;
    return 0;
  };

  double operator()( const vectord &x1, const vectord &x2)
  {
    double rl = computeScaledNorm2(x1,x2);
    double k = rl*rl;
    return exp(-k/2);
  };
  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    double rl = computeScaledNorm2(x1,x2);
    double k = rl*rl;
    return exp(-k/2)*k;
  };
};


////////////////////////////////////////////////////////////////////////////////
/** \brief Square exponential (Gaussian) kernel. ARD version. */
class SEArd: public ARDkernel
{
public:
  int init(size_t input_dim)
  {
    n_params = input_dim;
    n_inputs = input_dim;
    return 0;
  };

  double operator()( const vectord &x1, const vectord &x2 )
  {
    vectord ri = computeScaledDiff(x1,x2);
    double rl = norm_2(ri);
    double k = rl*rl;
    return exp(-k/2);
  };
  
  double gradient(const vectord &x1, const vectord &x2,
		  size_t component)
  {
    vectord ri = computeScaledDiff(x1,x2);
    double rl = norm_2(ri);
    double k = rl*rl;
    return exp(-k/2)*sqrt(ri(component));
  };
};


//@}

#endif
