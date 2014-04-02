/**  \file empiricalbayes.hpp \brief Empirical Bayesian estimator */
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


#ifndef  _EMPIRICALBAYES_HPP_
#define  _EMPIRICALBAYES_HPP_

#include <boost/scoped_ptr.hpp>
#include "criteria_functors.hpp"
#include "inneroptimization.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

 /**
   * \brief Bayesian optimization using different non-parametric 
   * processes as distributions over surrogate functions. 
   */
  class BAYESOPT_API EmpiricalBayes
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.h)
     */
    EmpiricalBayes(size_t dim, bopt_params params);

    /** 
     * Default destructor
     */
    virtual ~EmpiricalBayes();

   
    void updateHyperParameters();

    void setSamples(const matrixd &x, const vectord &y);
    void setSample(const vectord &x, double y);
    void addSample(const vectord &x, double y);


    Criteria* getCriteria();
    NonParametricProcess* getSurrogateModel();
    Dataset* getData();

  protected:
    bopt_params mParameters;                       ///< Configuration parameters
    size_t mDims;                                      ///< Number of dimensions
    randEngine mEngine;                             ///< Random number generator
    boost::scoped_ptr<Criteria> mCrit;                   ///< Metacriteria model
    boost::scoped_ptr<NonParametricProcess> mGP; ///< Pointer to surrogate model
    Dataset mData;                  ///< Dataset (x-> inputs, y-> labels/output)
    MeanModel mMean;

  private:

    EmpiricalBayes();

    void setSurrogateModel();    
    void setCriteria();

    /** 
     * \brief Checks the parameters and setups the elements (criteria, 
     * surrogate, etc.)
     */
    void __init__();

    CriteriaFactory mCFactory;
    boost::scoped_ptr<NLOPT_Optimization> kOptimizer;
  };

  /**@}*/


  inline vectord EmpiricalBayes::getPointAtMinimum() 
  { return mData.getPointAtMinimum(); };
  
  inline double EmpiricalBayes::getValueAtMinimum()
  { return mData.getValueAtMinimum(); };

  inline  Criteria* EmpiricalBayes::getCriteria()
  { return mCrit.get(); };

  inline  NonParametricProcess* EmpiricalBayes::getSurrogateModel()
  { return mGP.get(); };

  inline bopt_params* EmpiricalBayes::getParameters() 
  {return &mParameters;};


} //namespace bayesopt


#endif
