
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
     * @return 0 if terminate successfully, any other value otherwise
     */
    int optimize(vectord &bestPoint);

    /** 
     * \brief Execute ONE step the optimization process of the
     * function defined in evaluateSample.  
     * @param ii iteration number.
     */  
    void stepOptimization(size_t ii);

    /** Initialize the optimization process.  */
    void initializeOptimization();

    /** Once the optimization has been perfomed, return the optimal point. */
    virtual vectord getFinalResult() = 0;

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

    void setSamples(const matrixd &x, const vectord &y);
    void setSample(const vectord &x, double y);
    void addSample(const vectord &x, double y);
    vectord getPointAtMinimum();
    double getValueAtMinimum();

    NonParametricProcess* getSurrogateModel();
    void setSurrogateModel();    
    void  setCriteria();
    bopt_params* getParameters();
    randEngine& getRandomNumberGenerator();

  protected:
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

    /** Sample a single point in the input space. Used for epsilon greedy exploration. */
    virtual vectord samplePoint() = 0;

    /** 
     * \brief Selects the next point to evaluate according to a certain
     * criteria or metacriteria
     * 
     * @return next point to evaluate
     */
    vectord nextPoint();  

    void fitSurrogateModel();

    void plotInitialPoints();
  protected:
    Dataset mData;                  ///< Dataset (x-> inputs, y-> labels/output)
    boost::scoped_ptr<NonParametricProcess> mGP; ///< Pointer to surrogate model
    boost::scoped_ptr<Criteria> mCrit;                   ///< Metacriteria model
    bopt_params mParameters;                       ///< Configuration parameters
    size_t mDims;                                      ///< Number of dimensions
    randEngine mEngine;                             ///< Random number generator
    MeanModel mMean;

  private:

    BayesOptBase();

    /** 
     * \brief Checks the parameters and setups the elements (criteria, 
     * surrogate, etc.)
     */
    void __init__();

    CriteriaFactory mCFactory;
    NLOPT_Optimization* kOptimizer;
  };

  /**@}*/

  inline void BayesOptBase::setSamples(const matrixd &x, const vectord &y)
  { 
    mData.setSamples(x,y);  
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd 
  }

  inline void BayesOptBase::setSample(const vectord &x, double y)
  { 
    matrixd xx(1,x.size());  vectord yy(1);
    row(xx,0) = x;           yy(0) = y;
    mData.setSamples(xx,yy);   
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd
  }

  inline void BayesOptBase::addSample(const vectord &x, double y)
  {  mData.addSample(x,y); mMean.addNewPoint(x);  };

  inline vectord BayesOptBase::getPointAtMinimum() 
  { return mData.getPointAtMinimum(); };
  
  inline double BayesOptBase::getValueAtMinimum()
  { return mData.getValueAtMinimum(); };

  inline double BayesOptBase::evaluateCriteria( const vectord &query )
  {
    if (checkReachability(query))  return (*mCrit)(query);
    else return 0.0;
  };

  inline  NonParametricProcess* BayesOptBase::getSurrogateModel()
  { return mGP.get(); };

  inline bopt_params* BayesOptBase::getParameters() 
  {return &mParameters;};

  inline randEngine& BayesOptBase::getRandomNumberGenerator() 
  {return mEngine;};

  inline void BayesOptBase::fitSurrogateModel()
  {
    vectord optimalTheta = mGP->getHyperParameters();

    FILE_LOG(logDEBUG) << "Initial kernel parameters: " << optimalTheta;
    kOptimizer->run(optimalTheta);
    mGP->setHyperParameters(optimalTheta);
    FILE_LOG(logDEBUG) << "Final kernel parameters: " << optimalTheta;	

    mGP->fitSurrogateModel(); 
  };


} //namespace bayesopt


#endif
