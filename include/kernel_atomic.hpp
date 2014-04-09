/** \file kernel_atomic.hpp \brief Atomic (simple) kernel functions */
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

#ifndef  _KERNEL_ATOMIC_HPP_
#define  _KERNEL_ATOMIC_HPP_

#include <valarray>
#include "kernel_functors.hpp"
#include "elementwise_ublas.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{


  /** \brief Abstract class for an atomic kernel */
  class AtomicKernel : public Kernel
  {
  public:
    virtual void init(size_t input_dim)
    { n_inputs = input_dim; };

    void setHyperParameters(const vectord &theta) 
    {
      if(theta.size() != n_params)
	{
	  FILE_LOG(logERROR) << "Wrong number of kernel hyperparameters"; 
	  throw std::invalid_argument("Wrong number of kernel hyperparameters");
	}
      params = theta; //TODO: To make enough space. Make it more efficient.
      std::transform(theta.begin(), theta.end(), params.begin(), (double (*)(double)) exp);
      //      params = exp(theta);
    };

    vectord getHyperParameters() 
    { 
      vectord theta(params.size());
      std::transform(params.begin(), params.end(), theta.begin(), (double (*)(double)) log);
      return theta;
      //  return log(params);
    };
    size_t nHyperParameters() {return n_params;};

    virtual ~AtomicKernel(){};

  protected:
    size_t n_params;
    vectord params;
  };

  //////////////////////////////////////////////////////////////////////////
  /** \brief Constant kernel.  
   * Combined with sum and product kernel can be used to scale a kernel
   * or add noise.
   */

  class ConstKernel: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    { return params(0); };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { return 0.0; };
  };


  /** \brief Linear kernel. */
  class LinKernel: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 0;  n_inputs = input_dim; };

    double operator()(const vectord &x1, const vectord &x2)
    {
      assert(x1.size() == x2.size());
      return boost::numeric::ublas::inner_prod(x1,x2); 
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false);  return 0.0;   };
  };

  /** \brief Linear kernel. */
  class LinKernelARD: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      assert(x1.size() == x2.size());
      vectord v1 = utils::ublas_elementwise_div(x1, params);
      vectord v2 = utils::ublas_elementwise_div(x2, params);
      return boost::numeric::ublas::inner_prod(v1,v2); 
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false); return 0.0;  };
  };

  /** Kernel for categorical data. It measures the hamming distance
   *  between vectors. 
   */
  class HammingKernel: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    size_t hammingDistance(const vectori& s1, const vectori& s2)
    {
      size_t hdist = 0;
      vectori::const_iterator i1;
      vectori::const_iterator i2;
      for( i1 = s1.begin(), i2 = s2.begin(); 
	   i1 < s1.end() && i2 < s2.end();
	   ++i1, ++i2 )
	{
	  hdist += (*i1 == *i2) ? 0 : 1;
	}
      return hdist;
    }

    double operator()(const vectord &x1, const vectord &x2)
    { 
      const size_t n = x1.size();
      const double coef = -params(0)/2.0;
      vectori s1(n);
      vectori s2(n);

      for(size_t i=0; i<n; ++i)
	{
	  // We add 0.5 to avoid floating point approximation errors
	  s1(i) = static_cast<int>(x1(i)+0.5);
	  s2(i) = static_cast<int>(x2(i)+0.5);
	}
      
      const double dist = static_cast<double>(hammingDistance(s1,s2));
      return std::exp(coef*dist*dist);
    };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { return 0.0; };
  };


  ///////////////////////////////////////////////////////////////////////////
  /** \brief Abstract class for isotropic kernel functors */
  class ISOkernel : public AtomicKernel
  {
  public:
    virtual ~ISOkernel(){};

  protected:
    inline double computeWeightedNorm2(const vectord &x1, const vectord &x2)
    {  
      assert(n_inputs == x1.size());
      assert(x1.size() == x2.size());
      return norm_2(x1-x2)/params(0); 
    };
  };



  /** \brief Abstract class for anisotropic kernel functors. 
   * Typically ARD (Automatic Relevance Determination)
   */
  class ARDkernel : public AtomicKernel
  {
  public:
    virtual ~ARDkernel(){};

  protected:
    inline double computeWeightedNorm2(const vectord &x1, const vectord &x2)
    {
      assert(n_inputs == x1.size());
      assert(x1.size() == x2.size());
      assert(x1.size() == params.size());

      vectord xd = x1-x2;
      vectord r = utils::ublas_elementwise_div(xd, params);
      return norm_2(r);
    };
  };



  /** \brief Matern isotropic kernel of 1st order */
  class MaternIso1: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1;  n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      double r = computeWeightedNorm2(x1,x2);
      return exp(-r);
    };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    {
      double r = computeWeightedNorm2(x1,x2);
      return r*exp(-r);
    };
  };


  /** \brief Matern ARD kernel of 1st order */
  class MaternARD1: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim; n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      double r = computeWeightedNorm2(x1,x2);
      return exp(-r);
    };

    //TODO: 
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false);  return 0.0;  };
  };


  /** \brief Matern kernel of 3rd order */
  class MaternIso3: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r)*er;
    };

    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return r*r*er; 
    };
  };

  /** \brief Matern ARD kernel of 3rd order */
  class MaternARD3: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r)*er;
    };

    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {
      assert(false); return 0.0;
    };
  };


  /** \brief Matern isotropic kernel of 5th order */
  class MaternIso5: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r*(1+r/3))*er;
    };
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {    
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return r*(1+r)/3*r*er; 
    };
  };


  /** \brief Matern ARD kernel of 5th order */
  class MaternARD5: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r*(1+r/3))*er;
    };

    //TODO:
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    { assert(false); return 0.0; };
  };

  /** Polynomial covariance function*/
  class Polynomial: public AtomicKernel
  {
  public:
    Polynomial(){ mExp = 1; };

    void init(size_t input_dim)
    { n_params = 2;  n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double xx = boost::numeric::ublas::inner_prod(x1,x2); 
      return params(0)*params(0) * pow((params(1)+xx),mExp);
    };

    //TODO:
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    { assert(false); return 0.0; };
  protected:
    size_t mExp;
  };

  class Polynomial2: public Polynomial { public: Polynomial2(){ mExp = 2;};};
  class Polynomial3: public Polynomial { public: Polynomial3(){ mExp = 3;};};
  class Polynomial4: public Polynomial { public: Polynomial4(){ mExp = 4;};};
  class Polynomial5: public Polynomial { public: Polynomial5(){ mExp = 5;};};
  class Polynomial6: public Polynomial { public: Polynomial6(){ mExp = 6;};};
  class Polynomial7: public Polynomial { public: Polynomial7(){ mExp = 7;};};

  /** \brief Square exponential (Gaussian) kernel. Isotropic version. */
  class SEIso: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl;
      return exp(-k/2);
    };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl;
      return exp(-k/2)*k;
    };
  };


  /** \brief Square exponential (Gaussian) kernel. ARD version. */
  class SEArd: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2 )
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl;
      return exp(-k/2);
    };
  
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl;
      double r = (x1(component) - x2(component))/params(component);
      return exp(-k/2)*r*r;
    };
  };


  /** \brief Rational quadratic (Student's t) kernel. Isotropic version. */
  class RQIso: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 2; n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl/(2*params(1));
      return pow(1+k,-params(1));
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false); return 0.0;  };
  };

  //@}

} //namespace bayesopt

#endif
