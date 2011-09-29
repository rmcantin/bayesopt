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


#ifndef  _GAUSSPROCESS_HPP_
#define  _GAUSSPROCESS_HPP_

// BOOST Libraries
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include "krigwpr.h"

// Default values

//#define FIBONACCI_SEED  123u
//#define MT_SEED         156u
#define KERNEL_THETA    0.21
#define KERNEL_P        1.6
#define PRIOR_ALPHA     1.0
#define PRIOR_BETA      0.1
#define PRIOR_DELTA_SQ  10.0
#define DEF_REGULARIZER 1e-4
#define MAX_ITERATIONS  300
#define MAX_DIM         20

	
using namespace boost::numeric::ublas;	

typedef std::vector<vector<double> >  vecOfvec;
typedef symmetric_matrix<double,lower,row_major,bounded_array<double,90000> > covMatrix;
 
/** \addtogroup BayesOptimization */
/*@{*/



class GaussianProcess
{
public:

  GaussianProcess();  
  GaussianProcess( double theta, double p,
		   double alpha, double beta, 
		   double delta, double noise);
  GaussianProcess( gp_params params );
  ~GaussianProcess();

  void setSamples(matrix<double> x, vector<double> y)
  {
    size_t nPoints = x.size1();
    
    for (size_t i=0; i<nPoints; ++i)
      mGPXX.push_back(row(x,i));

    mGPY = y;
  }

  void addSample(vector<double> x, double y)
  {
    mGPXX.push_back(x);
    mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
  }

  /** 
   * Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * @param yPred mean of the predicted Gaussian distribution
   * @param sPred std of the predicted Gaussian distribution
   * 
   * @return error code.
   */	
  int prediction(const vector<double> &query,
		 double& yPred, double& sPred);

  //  int marginalLikelihood(const vector<double> &query,
			 
			 

  int fitGP();
  int precomputeGPParams();
  int addNewPointToGP( const vector<double> &Xnew,
		       double Ynew);

  vector<double> getPointAtMinimum()
  { return mGPXX[mMinIndex]; }

  double getValueAtMinimum()
  { return mGPY(mMinIndex); }

  /*double getNormValue()
    { return mYNorm(mMinIndex); }*/

  double getTheta()
  { return mTheta; }

  void setTheta( double theta )
  { mTheta = theta; }

protected:
  double correlationFunction( const vector<double> &x1,
			      const vector<double> &x2 );

  int computeCorrMatrix();

  vector<double> computeCrossCorrelation(const vector<double> &query);

  double meanFunction( const vector<double> &x);

  /*  inline void normalizeData()
  {
    scalar_vector<double> MinYVec(mGPY.size(), mGPY(mMinIndex));
    mYNorm = (mGPY - MinYVec) * ( 1/(mGPY(mMaxIndex)-mGPY(mMinIndex)) );
    };*/

  inline void checkBoundsY( size_t i )
  {
    if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
    else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
  };

protected:
  double mTheta, mP;            // Kernel parameters
  const double mAlpha, mBeta;         // GP prior parameters (Inv-Gamma)
  const double mDelta2,mRegularizer;  // GP prior parameters (Normal)

  vecOfvec mGPXX;                     // TODO:Data inputs
  vector<double> mGPY;                // Data values
  // vector<double> mYNorm;              // Normalized data values
	
  double mMu, mSig;                   // GP posterior paramethers
  // Precomputed GP prediction operations
  
  bounded_matrix<double, 
		 MAX_ITERATIONS, 
		 MAX_ITERATIONS> mInvR; // Inverse Correlation matrix
  

  //covMatrix mInvR;

  vector<double> mUInvR;              
  double mUInvRUDelta;
  vector<double> mYUmu;

  size_t mMinIndex, mMaxIndex;

  int mVerbose;
};


/**@}*/
// end namespaces

#endif
