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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "randgen.hpp"

#include "elementwiseUblas.hpp"
#include "krigwpr.h"

// Default values

//#define FIBONACCI_SEED  123u
//#define MT_SEED         156u
#define MAX_ITERATIONS  300
#define MAX_DIM         20

// DIRECT default values
#define MAX_DIRECT_EVALUATIONS  1000
#define MAX_DIRECT_ITERATIONS   300

// Latin Hypercube Sampling (LHS) default values
#define N_LHS_EVALS_PER_DIM     10
#define MAX_LHS_EVALUATIONS     100
	
using namespace boost::numeric::ublas;	


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
class Krigging
{

 public:
  
  /**
   *  Default Constructor 
   */
  Krigging();

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
  Krigging( double theta, double p,
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
  Krigging( gp_params params,
	    size_t nIter, bool useCool = false); 
	
  /** 
   * Default destructor
   * 
   * @return 
   */
  virtual ~Krigging();

  /** 
   * Execute the optimization process of the function defined in evaluateSample.
   * We assume that the function is defined in the [0,1] hypercube, as a normalized
   * representation of the bound constrains.
   * 
   * @see scaleInput
   * @see evaluateSample
   *
   * @param bestPoint returns the optimum value in a ublas::vector defined in the 
   *                  hypercube [0,1], it might also be used as an initial point
   * @param mtRandom random engine from boost random library
   * 
   * @return 1 if terminate successfully, 0 otherwise
   */
  int optimize( vector<double> &bestPoint,
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
  int optimize( vector<double> &bestPoint,
		vector<double> &lowerBound,
		vector<double> &upperBound,
		randEngine& mtRandom);

  /** 
   * Function that returns the negative Expected Improvement (-EI) of a series of queries
   * in the hypercube [0,1] in order to choose the best point to try the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * 
   * @return negative Expected Improvement (-EI).
   */	
  double negativeExpectedImprovement( const vector<double> &query );

  /** 
   * Function that returns the Lower Confidence Bound (LCB) of a series of queries
   * in the hypercube [0,1] in order to choose the best point to try the next iteration.
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * 
   * @return Lower Confidence Bound (LCB).
   */	
  double lowerConfidenceBound(const vector<double> &query);

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
  virtual double evaluateSample( const vector<double> &query ) 
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
  virtual bool checkReachability( const vector<double> &query )
  { return true; };
  /** 
   * Sets the parameters for the LCB criterium
   * LCB = mean - beta * std
   * 
   * @beta value of std coefficient
   */
  void setLCBparams(double beta)
  { mLCBbeta = beta; }

 /** 
   * Sets the parameters for the EI criterium
   * EI = E[max(y-y*,0)^g]
   * 
   * @g exponent of the improvement function
   */
  void setEIparams(double g)
  { mEIg = g; }

  /** 
   * Defines which cretirium to use. So far, only Expected Improvement 
   * or Lower Confidence Bound are defined.
   * 
   * @param useEI if true, use EI, if false, use LCB
   */
  void setCriteria(bool useEI)
  { mUseEI = useEI; }
 
protected:

  int allocateMatrices(size_t nSamples, size_t nDims);

  int sampleInitialPoints( size_t nSamples, 
			   size_t nDims,
			   bool useLatinBox,
			   randEngine& mtRandom);

  int checkBoundsY( size_t i);

  int updateCoolingScheme(size_t nTotalIterations,
			  size_t nCurrentIteration);
  int fitGP();
 	
  int addNewPointToGP( const vector<double> &Xnew );
	
  int nextPoint( vector<double> &Xnext );
  int nextPointDIRECT(vector<double> &Xnext);
  int nextPointNLOPT(vector<double> &Xnext);

  double correlationFunction( const vector<double> &x1,
			      const vector<double> &x2 );

  inline double evaluateNormalizedSample( const vector<double> &query);
  inline void normalizeData();

  // Math functions
  unsigned int factorial(unsigned int no, unsigned int a = 1);
  double pdf(double x);
  double cdf(double x);

protected:

  // Member variables
  bool mUseEI, mUseNLOPT;
  double mLCBbeta;                   // LCB = mean - param * std

  const double mTheta, mP;            // Kernel parameters
  const double mAlpha, mBeta;         // GP prior parameters (Inv-Gamma)
  const double mDelta2,mRegularizer;  // GP prior parameters (Normal)
  size_t mMaxIterations;
  const size_t mMaxDim;// Maximum Krigging evaluations and dimensions

  const bool mUseCool;
  unsigned int mEIg;

  vector<double> mLowerBound;
  vector<double> mRangeBound;

  //bounded_matrix<double, MAX_ITERATIONS, MAX_DIM> mGPX;  // Data points
  //bounded_vector<double, MAX_ITERATIONS>          mGPY;  // Data values
  //bounded_vector<double, MAX_ITERATIONS>          mYNorm;// Normalized data values

  matrix<double> mGPX;  // Data points
  vector<double> mGPY;  // Data values
  vector<double> mYNorm;// Normalized data values
	
  double mMu, mSig;                   // GP posterior paramethers
	
  double mMinY, mMaxY;
  vector<double> mMinX;

  bounded_matrix<double, MAX_ITERATIONS, MAX_ITERATIONS> mInvR; // Inverse Correlation matrix	
  vector<double> mUInvR;              // Precomputed GP prediction operations
  double mUInvRUDelta;
  vector<double> mYUmu;

  int mVerbose;

};

/**@}*/
// end namespaces


#endif
