/** \file nonparametricprocess.hpp 
    \brief Nonparametric process abstract module */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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


#ifndef __NONPARAMETRICPROCESS_HPP__
#define __NONPARAMETRICPROCESS_HPP__

#include "parameters.h"
#include "kernel_functors.hpp"
#include "mean_functors.hpp"
#include "randgen.hpp"
#include "specialtypes.hpp"
#include "inneroptimization.hpp"	

#define USE_CHOL 1


/** \addtogroup NonParametricProcesses */
/**@{*/

/**
 * \brief Abstract class to implement non-parametric processes
 */
class NonParametricProcess: public InnerOptimization
{
public:
  NonParametricProcess(double noise = DEFAULT_NOISE);
  
  virtual ~NonParametricProcess();

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
  virtual int prediction(const vectord &query,
			 double& yPred, double& sPred) = 0;

  /** 
   * Computes the negative log likelihood of the data.
   * 
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  virtual double negativeLogLikelihood(size_t index = 1) = 0;

  /** 
   * CURRENTLY NOT USED:
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  virtual double negativeLogLikelihood(double& grad,
				       size_t index = 1)
  {return 0.0;};

 
  /** 
   * Expected Improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * @param g exponent (used for annealing)
   * 
   * @return negative value of the expected improvement
   */
  virtual double negativeExpectedImprovement(const vectord &query,
					     size_t g = 1) = 0;

  /** 
   * Lower confindence bound. Can be seen as the inverse of the Upper 
   * confidence bound
   *
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param beta std coefficient (used for annealing)
   * 
   * @return value of the lower confidence bound
   */
  virtual double lowerConfidenceBound(const vectord &query,
				      double beta = 1) = 0;

  /** 
   * Probability of improvement algorithm for minimization
   * 
   * @param yPred mean of the prediction
   * @param sPred std of the prediction
   * @param yMin  minimum value found
   * @param epsilon minimum improvement margin
   * 
   * @return negative value of the probability of improvement
   */
  virtual double negativeProbabilityOfImprovement(const vectord &query,
						  double epsilon = 0.1) = 0;

  /** 
   * Sample outcome acording to the marginal distribution at the query point.
   * 
   * @param query query point
   * @param eng boost.random engine
   * 
   * @return outcome
   */
  virtual double sample_query(const vectord& query, randEngine& eng) = 0;

		 		 
  /** 
   *  Computes the GP based on mGPXX
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the begining
   * 
   * @return error code
   */
  int fitGP();
  
  /** 
   *  Add new point efficiently using Matrix Decomposition Lemma
   *  for the inversion of the correlation matrix. Maybe it is faster
   *  to just construct and invert a new matrix each time.
   * 
   * @return error code
   */   
  int addNewPointToGP( const vectord &Xnew,
		       double Ynew);


  // Getters and setters
  inline void setSamples(const matrixd &x, const vectord &y)
  {
    mGPY = y;
    for (size_t i=0; i<x.size1(); ++i)
      {
	mGPXX.push_back(row(x,i));
	checkBoundsY(i);
      } 
    mMeanV = (*mMean)(mGPXX);
  };

  inline void addSample(const vectord &x, double y)
  {
    mGPXX.push_back(x);
    mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
    checkBoundsY(mGPY.size()-1);
    mMeanV.resize(mMeanV.size()+1);  
    mMeanV(mMeanV.size()-1) = mMean->getMean(x);
  };

  inline double getSample(size_t index, vectord &x)
  {
    x = mGPXX[index];
    return mGPY(index);
  }

  inline vectord getPointAtMinimum()
  { return mGPXX[mMinIndex]; };

  inline double getValueAtMinimum()
  { return mGPY(mMinIndex); };

  /** 
   * Select kernel (covariance function) for the surrogate process.
   * 
   * @param thetav kernel parameters
   * @param k_name kernel name
   * @return error_code
   */
  int setKernel (const vectord &thetav,
		 kernel_name k_name);

  /** 
   * Wrapper of setKernel for c arrays
   */
  int setKernel (const double *theta, size_t n_theta,
		 kernel_name k_name)
  {
    vectord th(n_theta);
    std::copy(theta, theta+n_theta, th.begin());
    return setKernel(th, k_name);
  };

  /** 
   * Select the parametric part of the surrogate process.
   * 
   * @param muv mean function parameters
   * @param m_name mean function name
   * @return error_code
   */
  int setMean (const vectord &muv,
	       mean_name m_name);

  /** 
   * Wrapper of setMean for c arrays
   */
  int setMean (const double *mu, size_t n_mu,
	       mean_name m_name)
  {
    vectord vmu(n_mu);
    std::copy(mu, mu+n_mu, vmu.begin());
    return setMean(vmu, m_name);
  };


protected:

  double innerEvaluate(const vectord& query)
  { 
    mKernel->setScale(query);
    return negativeLogLikelihood(1);
  };

  //inline double meanFunction( const vectord &x, size_t param_index = 0 )  
  //{ return (*mMean)(x); };

  /** 
   * Precompute some values of the prediction that do not depends on the query
   * 
   * @return error code
   */
  virtual int precomputePrediction()
  {return 1;}


  /** 
   * Computes the inverse of the Correlation (Kernel or Gram) matrix  
   * 
   * @return error code
   */
  int computeInverseCorrelation();
  int addNewPointToInverse(const vectord& correlation,
			   double selfcorrelation);

  int computeCholeskyCorrelation();
  int addNewPointToCholesky(const vectord& correlation,
			    double selfcorrelation);

  matrixd computeCorrMatrix(int dth_index = -1);
  vectord computeCrossCorrelation(const vectord &query);


  inline void checkBoundsY( size_t i )
  {
    if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
    else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
  };


protected:
  const double mRegularizer;  ///< GP likelihood parameters (Normal)
  vecOfvec mGPXX;                                   ///< Data inputs
  vectord mGPY;                                     ///< Data values

  Kernel* mKernel;                   ///< Pointer to kernel function
  ParametricFunction* mMean;           ///< Pointer to mean function
  bool kernelSet, meanSet;

  size_t mMinIndex, mMaxIndex;	

  // Precomputed GP prediction operations
  covMatrix mInvR;                   ///< Inverse Correlation matrix



  //  vectord mTheta;                             ///< Kernel parameters


  matrixd mL;
  vectord mAlphaV;
  vectord mMeanV;              

};

/**@}*/

#endif
