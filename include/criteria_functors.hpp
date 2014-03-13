/**  \file criteria_functors.hpp \brief Criteria functions */
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

#ifndef  _CRITERIA_FUNCTORS_HPP_
#define  _CRITERIA_FUNCTORS_HPP_

#include <map>
#include <algorithm>
#include "optimizable.hpp"
#include "bayesianregressor.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions
   * \brief Set of criterium functions to select the next point during
   * optimization/exploration
   */
  //@{

  /**
   * \brief Abstract interface for criteria functors.
   */
  class Criteria: public RBOptimizable
  {
  public:
    virtual ~Criteria() {};
    virtual int init(BayesianRegressor *proc) { return 0; };
    virtual int init(BayesianRegressor *proc, 
		     const std::vector<Criteria*>& list) { return 0; };

    virtual bool requireComparison() = 0;
    double evaluate(const vectord &x) {return (*this)(x);}
    virtual double operator()(const vectord &x) = 0;
    virtual std::string name() = 0;
    virtual int setParameters(const vectord &params) = 0;
    virtual size_t nParameters() = 0;

    //Dummy functions.
    virtual void reset() { assert(false); };
    virtual bool checkIfBest(vectord& xNext,
			     std::string& name,
			     int& error_code) { assert(false); return false; };
  protected:
    BayesianRegressor *mProc;
  };


  template <typename CriteriaType> Criteria * create_func()
  {
    return new CriteriaType();
  }


  /** 
   * \brief Factory model for criterion functions
   * This factory is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   */
  class CriteriaFactory
  {
  public:
    CriteriaFactory ();
    virtual ~CriteriaFactory () {};
  
    //Criteria* create(criterium_name name, BayesianRegressor* proc);
    Criteria* create(std::string name, BayesianRegressor* proc);
    
  private:
    typedef Criteria* (*create_func_definition)();
    std::map<std::string , CriteriaFactory::create_func_definition> registry;
  };


  //@}

} //namespace bayesopt


#endif
