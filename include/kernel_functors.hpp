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

#include <map>
#include "parameters.h"
#include "specialtypes.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions
   * \brief Set of kernel or covariance functions for the nonparametric
   * processes.
   */
  //@{

  /** \brief Interface for kernel functors */
  class Kernel
  {
  public:
    virtual int init(size_t input_dim) {return 0;};
    virtual int init(size_t input_dim, Kernel* left, Kernel* right) {return 0;};

    virtual int setHyperParameters(const vectord &theta) = 0;
    virtual vectord getHyperParameters() = 0;
    virtual size_t nHyperParameters() = 0;

    virtual double operator()( const vectord &x1, const vectord &x2 ) = 0;
    virtual double gradient( const vectord &x1, const vectord &x2,
			     size_t component ) = 0;
    virtual ~Kernel(){};

  protected:
    size_t n_inputs;
  };

  template <typename KernelType> Kernel * create_func()
  {
    return new KernelType();
  }


  /** 
   * \brief Factory model for kernel functions
   * This factory is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   */
  class KernelFactory
  {
  public:
    KernelFactory ();
    virtual ~KernelFactory () {};
  
    Kernel* create(kernel_name name, size_t input_dim);
    Kernel* create(std::string name, size_t input_dim);
    
  private:
    typedef Kernel* (*create_func_definition)();
    std::map<std::string , KernelFactory::create_func_definition> registry;
  };

  //@}

} //namespace bayesopt


#endif
