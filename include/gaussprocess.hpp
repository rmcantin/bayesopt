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

#include "krigwpr.h"
#include "specialtypes.hpp"
	
//using namespace boost::numeric::ublas;	

 
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

  void setSamples(matrixd x, vectord y)
  {
    size_t nPoints = x.size1();
    
    for (size_t i=0; i<nPoints; ++i)
      mGPXX.push_back(row(x,i));

    mGPY = y;
  }

  void addSample(vectord x, double y)
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
  int prediction(const vectord &query,
		 double& yPred, double& sPred);

  //  int marginalLikelihood(const vectord &query,
			 
			 

  int fitGP();
  int precomputeGPParams();
  int addNewPointToGP( const vectord &Xnew,
		       double Ynew);

  vectord getPointAtMinimum()
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
  double correlationFunction( const vectord &x1,
			      const vectord &x2 );


  double correlationFunction( const vectord &x1, 
			      const vectord &x2,
			      double param, double &grad);


  int computeCorrMatrix();

  vectord computeCrossCorrelation(const vectord &query);

  double meanFunction( const vectord &x);

  /*  inline void normalizeData()
  {
    svectord MinYVec(mGPY.size(), mGPY(mMinIndex));
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
  vectord mGPY;                // Data values
  // vectord mYNorm;              // Normalized data values
	
  double mMu, mSig;                   // GP posterior paramethers
  


  // Precomputed GP prediction operations
  covMatrix mInvR;                   // Inverse Correlation matrix

  vectord mUInvR;              
  double mUInvRUDelta;
  vectord mYUmu;

  size_t mMinIndex, mMaxIndex;

  int mVerbose;
};


/**@}*/
// end namespaces

#endif
