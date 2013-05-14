/**  \file bayesoptbase.hpp \brief Bayesian optimization module */
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


#ifndef  _BAYESOPTBASE_HPP_
#define  _BAYESOPTBASE_HPP_

#include "criteria_functors.hpp"

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
  class BAYESOPT_API BayesOptBase
  {
  public:
  
    /** 
     * Default constructor
     */
    BayesOptBase();
    
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
     * \brief Execute the optimization process of the function defined
     * in evaluateSample.
     * 
     * @see evaluateSample
     * @see checkReachability
     *
     * @param bestPoint returns point with the optimum value in a ublas::vector.
     * @return 0 if terminate successfully, any other value otherwise
     */
    int optimize(vectord &bestPoint);

    /** 
     * \brief Execute ONE step the optimization process of the function defined
     * in evaluateSample.
     * 
     * @see evaluateSample
     * @see checkReachability
     *
     * @param ii iteration number.
     * @return 0 if terminate successfully, any other value otherwise
     */  
    int stepOptimization(size_t ii);

    /** 
     * Initialize the optimization process.
     * @return error_code
     */
    virtual int initializeOptimization() = 0;

    /** 
     * Once the optimization has been perfomed, return the optimal point.
     */
    virtual vectord getFinalResult() = 0;


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
     * \brief Evaluate the criteria considering if the query is
     * reachable or not.  This is a way to include non-linear
     * restrictions.
     *
     * @see checkReachability
     *
     * @param query query point
     * @return value of the criteria, 0 otherwise.
     */
    inline double evaluateCriteria( const vectord &query )
    {
      bool reachable = checkReachability(query);
      if (!reachable)  return 0.0;
      return (*mCrit)(query);
    };

    
    inline NonParametricProcess* getSurrogateModel()
    { return mGP.get(); };

    int setSurrogateModel();    
    int setCriteria();

  protected:
    /** 
     * Print data for every step according to the verbose level
     * 
     * @param iteration 
     * @param xNext 
     * @param yNext 
     * 
     * @return error code
     */
    virtual int plotStepData(size_t iteration, const vectord& xNext,
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
     * @return error code
     */
    virtual int findOptimal(vectord &xOpt) = 0;
  
    /** 
     * \brief Selects the initial set of points to build the surrogate
     * model.
     * @return error code
     */
    virtual int sampleInitialPoints() = 0;

    /** 
     * \brief Selects the next point to evaluate according to a certain
     * criteria or metacriteria
     * 
     * @param Xnext next point to evaluate
     * @return error code
     */
    int nextPoint( vectord &Xnext );  



  protected:
    boost::scoped_ptr<NonParametricProcess> mGP;  ///< Pointer to surrogate model
    boost::scoped_ptr<Criteria> mCrit;                    ///< Metacriteria model
    bopt_params mParameters;                        ///< Configuration parameters
    size_t mDims;                                       ///< Number of dimensions

  private:

    /** 
     * \brief Checks the parameters and setups the elements (criteria, 
     * surrogate, etc.)
     */
    int __init__();

    CriteriaFactory mCFactory;
    randEngine mEngine;
  };

  /**@}*/

} //namespace bayesopt


#endif
