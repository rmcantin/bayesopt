/**  \file bayesopt.hpp \brief Bayesian optimization C++-API*/
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

#ifndef  _BAYESOPTAPI_HPP_
#define  _BAYESOPTAPI_HPP_

#include "boundingbox.hpp"
#include "bayesoptbase.hpp"
#include "inneroptimization.hpp"

namespace bayesopt  {

  /** \addtogroup BayesOpt */
  /**@{*/

  /**
   * \brief Bayesian optimization using different non-parametric
   * processes as distributions over surrogate functions. The
   * exploration spaces is assumed to be continous and box-bounded.
   */
  class BAYESOPT_API ContinuousModel: public BayesOptBase
  {
  public:
   
    /** Default constructor */
    //    ContinuousModel();

    /** 
     * Constructor
     * @param dim number of input dimensions
     * @param params set of parameters (see parameters.h)
     */
    ContinuousModel(size_t dim, bopt_params params);

    /**  Default destructor  */
    virtual ~ContinuousModel();
  
    /** 
     * Once the optimization has been perfomed, return the optimal
     * point.
     */
    vectord getFinalResult();

    /** 
     * \brief Sets the bounding box. 
     *
     * @param lowerBound vector with the lower bounds of the hypercube
     * @param upperBound vector with the upper bounds of the hypercube
     */
    void setBoundingBox( const vectord &lowerBound,
			const vectord &upperBound);


  protected:

    /** 
     * \brief Print data for every step according to the verbose level
     * 
     * @param iteration iteration number 
     * @param xNext next point
     * @param yNext function value at next point
     */
    void plotStepData(size_t iteration, const vectord& xNext,
		      double yNext);

    /** Selects the initial set of points to build the surrogate model. */
    void sampleInitialPoints(matrixd& xPoints, vectord& yPoints);

    /** Sample a single point in the input space. Used for epsilon greedy exploration. */
    vectord samplePoint();

    /** 
     * \brief Wrapper for the target function normalize in the hypercube
     * [0,1]
     * @param query point to evaluate in [0,1] hypercube
     * @return actual return value of the target function
     */
    double evaluateSampleInternal( const vectord &query );
    
    /** 
     * \brief Wrapper of the innerOptimization class to find the optimal
     * point acording to the criteria.
     * @param xOpt optimal point
     */
    void findOptimal(vectord &xOpt);

  private:
    boost::scoped_ptr<utils::BoundingBox<vectord> > mBB;      ///< Bounding Box (input space limits)
    NLOPT_Optimization* cOptimizer;

    ContinuousModel();                       ///< Default constructor forbidden.
  };
  

  /**
   * \brief Sequential Kriging Optimization using different non-parametric 
   * processes as surrogate (kriging) functions. 
   */
  class BAYESOPT_API DiscreteModel : public BayesOptBase
  {
  public:
    /** 
     * Constructor
     * @param validSet  Set of potential inputs
     * @param params set of parameters (see parameters.h)
     */
    DiscreteModel(const vecOfvec &validSet, bopt_params params);
    
    /** Default destructor  */
    virtual ~DiscreteModel();

    /** Once the optimization has been perfomed, return the optimal point. */
    vectord getFinalResult();

    
  protected:
    
    
    /** Print data for every step according to the verbose level */
    void plotStepData(size_t iteration, const vectord& xNext,
		     double yNext);

    /** Selects the initial set of points to build the surrogate model. */
    void sampleInitialPoints(matrixd& xPoints, vectord& yPoints);

    /** Sample a single point in the input space. Used for epsilon greedy exploration. */
    vectord samplePoint();

    /** 
     * \brief Wrapper for the target function normalize in the hypercube
     * [0,1]
     * @param query point to evaluate in [0,1] hypercube
     * @return actual return value of the target function
     */
    double evaluateSampleInternal( const vectord &query ); 

    void findOptimal(vectord &xOpt);

  protected:
    vecOfvec mInputSet;               ///< List of input points

  private:
    DiscreteModel();         ///< Default constructor forbidden.
  };


  /**@}*/

  //////////////////////////////////////////////////////////////////////
  //                         Inline methods 
  //////////////////////////////////////////////////////////////////////

  inline double ContinuousModel::evaluateSampleInternal( const vectord &query )
  { return evaluateSample(mBB->unnormalizeVector(query));  };

  inline void ContinuousModel::findOptimal(vectord &xOpt)
  { cOptimizer->run(xOpt); };

  inline vectord ContinuousModel::samplePoint()
  {	    
    randFloat drawSample(mEngine,realUniformDist(0,1));
    vectord Xnext(mDims);    
    for(vectord::iterator x = Xnext.begin(); x != Xnext.end(); ++x)
      {	*x = drawSample(); }
    return Xnext;
  };



  inline vectord DiscreteModel::samplePoint()
  {   
    randInt sample(mEngine, intUniformDist(0,mInputSet.size()-1));
    return mInputSet[sample()];
  };

  inline double DiscreteModel::evaluateSampleInternal( const vectord &query )
  { return evaluateSample(query); }; 




}  //namespace bayesopt


#endif
