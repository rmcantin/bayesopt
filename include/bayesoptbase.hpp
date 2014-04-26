
/**  \file bayesoptbase.hpp \brief BayesOpt common module for interfaces */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2014 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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


#ifndef  _BAYESOPTBASE_HPP_
#define  _BAYESOPTBASE_HPP_

#include <boost/scoped_ptr.hpp>
#include <boost/random.hpp>
#include "parameters.h"
#include "specialtypes.hpp"
//#include "posteriormodel.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {


  //Forward declaration
  class PosteriorModel;
  class ProbabilityDistribution;
  class Dataset;

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

  /**
   * \brief Abstract module for Bayesian optimization.
   *
   * This module provides Bayesian optimization using different
   * non-parametric processes (Gaussian process or Student's t
   * process) as distributions over surrogate functions.
   *
   * \see ContinuousModel for implementations of this module for
   * a continuous input spaces
   *
   * \see DiscreteModel for implementations of this module for
   * a discrete input spaces or categorical input variables
   */
  class BAYESOPT_API BayesOptBase
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.h)
     */
    BayesOptBase(size_t dim, bopt_params params);

    /** 
     * Default destructor
     */
    virtual ~BayesOptBase();

    /** 
     * \brief Function that defines the actual function to be optimized.
     * This function must be modified (overriden) according to the
     * specific problem.
     *
     * @param query point to be evaluated. 
     * @return value of the function at the point evaluated.
     */
    virtual double evaluateSample( const vectord &query ) = 0;
    

    /** 
     * \brief This function checks if the query is valid or not. It can
     * be used to introduce arbitrary constrains. Since the Gaussian
     * process assumes smoothness, constrains are managed by the inner
     * optimizer (e.g.:DIRECT), being highly time consuming. If the
     * constrain is very tricky, DIRECT will need much more function
     * evaluations.
     *
     * Note: This function is experimental. Thus it is not made pure virtual.
     * Using it is completely optional.
     * 
     * @param query point to be evaluated.
     * 
     * @return boolean value showing if the the function is valid at
     *         the query point or not.
     */ 
    virtual bool checkReachability( const vectord &query )
    { return true; };


    /** 
     * \brief Execute the optimization process of the function defined
     * in evaluateSample.
     * 
     * @see evaluateSample
     * @see checkReachability
     *
     * @param bestPoint returns point with the optimum value in a ublas::vector.
     */
    void optimize(vectord &bestPoint);

    /** 
     * \brief Execute ONE step the optimization process of the
     * function defined in evaluateSample.  
     */  
    void stepOptimization();

    /** Initialize the optimization process.  */
    void initializeOptimization();

    /** Once the optimization has been perfomed, return the optimal point. */
    virtual vectord getFinalResult() = 0;

    ProbabilityDistribution* getPrediction(const vectord& query);
    const Dataset* getData();
    bopt_params* getParameters();
    double getValueAtMinimum();
    double evaluateCriteria(const vectord& query);

  protected:
    vectord getPointAtMinimum();

    /** 
     * Print data for every step according to the verbose level
     * 
     * @param iteration 
     * @param xNext 
     * @param yNext 
     */
    virtual void plotStepData(size_t iteration, const vectord& xNext,
			      double yNext) = 0;


    /** 
     * \brief Wrapper for the target function normalize in the hypercube
     * [0,1]
     * @param query point to evaluate in [0,1] hypercube
     * @return actual return value of the target function
     */
    virtual double evaluateSampleInternal( const vectord &query ) = 0;


    /** 
     * \brief Returns the optimal point acording to certain criteria
     * @see evaluateCriteria
     *
     * @param xOpt optimal point
     */
    virtual void findOptimal(vectord &xOpt) = 0;
  
    /** Selects the initial set of points to build the surrogate model. */
    virtual void sampleInitialPoints(matrixd& xPoints, vectord& yPoints) = 0;

    /** Sample a single point in the input space. Used for epsilon
	greedy exploration. */
    virtual vectord samplePoint() = 0;

  protected:
    bopt_params mParameters;                    ///< Configuration parameters
    size_t mDims;                                   ///< Number of dimensions
    size_t mCurrentIter;                        ///< Current iteration number
    boost::mt19937 mEngine;                      ///< Random number generator

  private:
    boost::scoped_ptr<PosteriorModel> mModel;

  private:

    BayesOptBase();

    /** 
     * \brief Selects the next point to evaluate according to a certain
     * criteria or metacriteria
     * 
     * @return next point to evaluate
     */
    vectord nextPoint();  

  };

  /**@}*/



} //namespace bayesopt


#endif
