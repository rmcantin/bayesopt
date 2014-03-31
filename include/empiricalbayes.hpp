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
#include "log.hpp"


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
  class EmpiricalBayes
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

    /** 
     * \brief Evaluate the criteria considering if the query is
     * reachable or not.  This is a way to include non-linear
     * restrictions.
     *
     * @see checkReachability
     *
     * @param query query point
     * @return value of the criteria, 0 otherwise.
     */
    double evaluateCriteria( const vectord &query );

    NonParametricProcess* getSurrogateModel();
    void setSurrogateModel();    
    void  setCriteria();

  protected:
 
    void fitSurrogateModel();


  protected:
    size_t mDims;                                      ///< Number of dimensions
    bopt_params mParameters;                       ///< Configuration parameters
    randEngine& mEngine;                            ///< Random number generator
    Dataset mData;                  ///< Dataset (x-> inputs, y-> labels/output)
    MeanModel mMean;
    boost::scoped_ptr<NonParametricProcess> mGP; ///< Pointer to surrogate model
    boost::scoped_ptr<Criteria> mCrit;                   ///< Metacriteria model

  private:
    CriteriaFactory mCFactory;
    NLOPT_Optimization* kOptimizer;
  };

  /**@}*/

  inline void EmpiricalBayes::setSamples(const matrixd &x, const vectord &y)
  { mData.setSamples(x,y);  mGP->setSamples(x,y);  }

  inline void EmpiricalBayes::setSample(const vectord &x, double y)
  { 
    matrixd xx(1,x.size());  vectord yy(1);
    row(xx,0) = x;           yy(0) = y;
    mData.setSamples(xx,yy);   mGP->setSamples(xx,yy);  }

  inline void EmpiricalBayes::addSample(const vectord &x, double y)
  {  mData.addSample(x,y); mGP->addSample(x,y);  };

  inline vectord EmpiricalBayes::getPointAtMinimum() 
  { return mData.getPointAtMinimum(); };
  
  inline double EmpiricalBayes::getValueAtMinimum()
  { return mData.getValueAtMinimum(); };

  inline double EmpiricalBayes::evaluateCriteria( const vectord &query )
  {
    if (checkReachability(query))  return (*mCrit)(query);
    else return 0.0;
  };

  inline NonParametricProcess* EmpiricalBayes::getSurrogateModel()
  { return mGP.get(); };

  inline bopt_params* EmpiricalBayes::getParameters() 
  {return &mParameters;};

  inline randEngine& EmpiricalBayes::getRandomNumberGenerator() 
  {return mEngine;};

  inline void EmpiricalBayes::fitSurrogateModel()
  {
    vectord optimalTheta = mGP->getHyperParameters();

    FILE_LOG(logDEBUG) << "Initial kernel parameters: " << optimalTheta;
    kOptimizer->run(optimalTheta);
    mGP->setHyperParameters(optimalTheta);
    FILE_LOG(logDEBUG) << "Final kernel parameters: " << optimalTheta;	

    mGP->fitSurrogateModel(); 
  };


} //namespace bayesopt
