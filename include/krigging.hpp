/**
 * @file   krigging.hpp
 * @author Ruben Martinez-Cantin <rmcantin@ist.isr.utl.pt>
 * @date   Thu Mar 26 02:12:36 2009
 * 
 * @brief  Efficient Global Optimization with hyperpriors.
 *
 * This is an efficient, C++ implementation of the Bayesian optimization
 * algorithm presented in the papers:
 *
 * ----
 * Ruben Martinez-Cantin, Nando de Freitas, Arnaud Doucet and Jose Castellanos.
 * Active Policy Learning for Robot Planning and Exploration under Uncertainty. 
 * Robotics: Science and Systems. 2007
 *
 * Ruben Martinez-Cantin, Nando de Freitas, Eric Brochu, Jose Castellanos and 
 * Arnaud Doucet (2009) A Bayesian Exploration-Exploitation Approach for Optimal
 * Online Sensing and Planning with a Visually Guided Mobile Robot. Autonomous 
 * Robots - Special Issue on Robot Learning, Part B, 27(3):93-103.
 * ----
 * 
 * Basically, it uses the active learning strategy to optimize an "arbitrary" 
 * funtion using few iterations.
 * 
 * Copyright: See COPYING file that comes with this distribution
 */


#ifndef  _KRIGGING_HPP_
#define  _KRIGGING_HPP_

// BOOST Libraries
#include "specialtypes.hpp"
#include "randgen.hpp"
#include "elementwiseUblas.hpp"

#include "inneroptimization.hpp"
#include "basicgaussprocess.hpp"

#include "krigwpr.h"
#include "criteria.hpp"

/** \addtogroup BayesOptimization */
/*@{*/



/** 
 * @brief Object that integrates a Bayesian optimization algorithm. 
 * This is an efficient, C++ implementation of the Bayesian optimization
 * algorithm presented in the papers:
 *
 * ----
 * Ruben Martinez-Cantin, Nando de Freitas, Arnaud Doucet and Jose Castellanos.
 * Active Policy Learning for Robot Planning and Exploration under Uncertainty. 
 * Robotics: Science and Systems. 2007
 *
 * Ruben Martinez-Cantin, Nando de Freitas, Eric Brochu, Jose Castellanos and 
 * Arnaud Doucet (2009) A Bayesian Exploration-Exploitation Approach for Optimal
 * Online Sensing and Planning with a Visually Guided Mobile Robot. Autonomous 
 * Robots - Special Issue on Robot Learning, Part B, 27(3):93-103.
 * ----
 * 
 * Basically, it uses the active learning strategy to optimize an "arbitrary" 
 * funtion using few iterations.
 * 
 */
class SKO: public InnerOptimization
{

 public:
  
  /**
   *  Default Constructor 
   */
  SKO();

  /** 
   * Constructor
   * 
   * @param theta        kernel bandwidth
   * @param p            kernel exponent (not used)
   * @param alpha        inverse gamma prior
   * @param beta         inverse gamma prior
   * @param delta        normal prior
   * @param noise        observation noise
   * @param nIter        number of iterations before stopping 
   * @param useCool      select Sasena cooling/annealing strategy
   */
  SKO( double theta, double p,
       double alpha, double beta, 
       double delta, double noise,
       size_t nIter, bool useCool = false); 

  /** 
   * Constructor
   * 
   * @param params structure with the GP parameters
   *     theta        kernel bandwidth
   *     p            kernel exponent (not used)
   *     alpha        inverse gamma prior
   *     beta         inverse gamma prior
   *     delta        normal prior
   *     noise        observation noise
   * @param nIter        number of iterations before stopping 
   * @param useCool      select Sasena cooling/annealing strategy
   */
  SKO( gp_params params,
       size_t nIter, bool useCool = false); 
	
  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~SKO();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the [0,1] hypercube, as a 
   * normalized representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector defined in 
   * the hypercube [0,1], it might also be used as an initial point
   * @param mtRandom random engine from boost random library
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  int optimize( vectord &bestPoint,
		randEngine& mtRandom);


  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the hypercube defined by 
   * the lower and upper vectors. 
   *
   * @param bestPoint returns the optimum value in a ublas::vector x, 
   * it might also be used as an initial point
   * @param lowerBound vector with the lower bounds of the hypercube 
   * @param upperBound vector with the upper bounds of the hypercube 
   * @param mtRandom random engine from boost random library
   * @param useEI to decide whether we use the EI or UCB criterium
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  int optimize( vectord &bestPoint,
		vectord &lowerBound,
		vectord &upperBound,
		randEngine& mtRandom);

  /** 
   * Function that defines the actual mathematical function to be optimized.
   *
   * Virtual function for polymorphism. 
   *
   * This function must need to be modified according to the specific problem.
   *
   * @param query point to be evaluated. It is automatically generated using the
   *              Expected Improvement algorithm.
   * 
   * @return value of the function at the point evaluated
   */
  virtual double evaluateSample( const vectord &query ) 
  { return 0.0; };

  /** 
   * This function checks if the query is valid or not. It can be used 
   * to introduce arbitrary constrains. Since the Gaussian process 
   * assumes smoothness, constrains are managed by DIRECT, being highly
   * time consuming. If the constrain is very tricky, DIRECT will need
   * much more function evaluations.
   *
   * Note: This function is experimental. 
   * 
   * @param query point to be evaluated.
   * 
   * @return boolean value showing if the the function is valid at
   *         the query point or not.
   */ 
  virtual bool checkReachability( const vectord &query )
  { return true; };

  /** 
   * Function that returns the corresponding criteria of a series 
   * of queries in the hypercube [0,1] in order to choose the best point to try
   * the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * 
   * @return negative criteria (Expected Improvement, LCB, A-optimality, etc.).
   */	

  inline double innerEvaluate( const vectord &query, vectord &grad )
  {return evaluateCriteria(query);}

protected:


  double evaluateCriteria( const vectord &query );


  inline int nextPoint(vectord &Xnext)
  {return innerOptimize(Xnext);}


  int allocateMatrices(size_t nSamples, size_t nDims);

  int sampleInitialPoints( size_t nSamples, 
			   size_t nDims,
			   bool useLatinBox,
			   randEngine& mtRandom);

  inline double evaluateNormalizedSample( const vectord &query)
  { 
    vectord unnormalizedQuery = ublas_elementwise_prod(query,
							    mRangeBound);
  
    unnormalizedQuery = ublas_elementwise_add(unnormalizedQuery,
					      mLowerBound);
    
    return evaluateSample(unnormalizedQuery);
  } // evaluateNormalizedSample

protected:

  // Member variables
  BasicGaussianProcess mGP;
  Criteria crit;

  size_t mMaxIterations;
  const size_t mMaxDim;// Maximum SKO evaluations and dimensions

  vectord mLowerBound;
  vectord mRangeBound;
	
  int mVerbose;

};

/**@}*/
// end namespaces


#endif
