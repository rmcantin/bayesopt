/** \file nonparametricprocess.hpp 
    \brief Nonparametric process abstract module */
/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#ifndef __NONPARAMETRICPROCESS_HPP__
#define __NONPARAMETRICPROCESS_HPP__

#include "parameters.h"
#include "kernel_functors.hpp"
#include "meanfuncs.hpp"
#include "randgen.hpp"
#include "specialtypes.hpp"
#include "inneroptimization.hpp"	


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
   * @param grad gradient of the negative Log Likelihood
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
  virtual double negativeExpectedImprovement(double yPred, double sPred,
					     double yMin, size_t g = 1)
  {return 0.0;};

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
  virtual double lowerConfidenceBound(double yPred, double sPred,
				     double beta = 1)
  {return 0.0;};

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
  virtual double negativeProbabilityOfImprovement(double yPred, double sPred,
						  double yMin, 
						  double epsilon = 0.1)
  {return 0.0;};
		 		 
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
  inline void setSamples(matrixd x, vectord y)
  {
    for (size_t i=0; i<x.size1(); ++i)
      mGPXX.push_back(row(x,i));

    mGPY = y;
  };

  inline void addSample(vectord x, double y)
  {
    mGPXX.push_back(x);
    mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
  };

  inline vectord getPointAtMinimum()
  { return mGPXX[mMinIndex]; };

  inline double getValueAtMinimum()
  { return mGPY(mMinIndex); };

  /*
  inline vectord getTheta()
  { return mTheta; };
  
  inline void setTheta( const vectord& theta )
  { mTheta = theta; };
  */

  /** 
   * Chooses which kernel to use in the surrogate process.
   * 
   * @param theta kernel parameter
   * @param k_name kernel name
   */
  int setKernel (double theta,
		 kernel_name k_name)
  {
    return setKernel(svectord(1,theta), k_name);
  };

  int setKernel (const vectord &thetav,
		 kernel_name k_name);

  virtual double sample_query(const vectord& query, 
			      randEngine& eng)
  { return 0.0; }

protected:
  double innerEvaluate(const vectord& query)
  { 
    mKernel->setScale(query);
    return negativeLogLikelihood(1);
  };

  inline double meanFunction( const vectord &x, size_t param_index = 0 )  
  { //TODO: Parametrize this
    return means::One(x); };

  /** 
   * Precompute some values of the prediction that do not depends on the query
   * 
   * @return error code
   */
  virtual int precomputeGPParams()
  {return 1;}


  /** 
   * Computes the inverse of the Correlation (Kernel or Gram) matrix  
   * 
   * @return error code
   */
  int computeInverseCorrMatrix();

  matrixd computeCorrMatrix(int dth_index = -1);

  vectord computeCrossCorrelation(const vectord &query);


  inline void checkBoundsY( size_t i )
  {
    if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
    else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
  };


protected:
  vecOfvec mGPXX;                                   ///< Data inputs
  vectord mGPY;                                     ///< Data values
	
  // Precomputed GP prediction operations
  covMatrix mInvR;                   ///< Inverse Correlation matrix

  size_t mMinIndex, mMaxIndex;

  Kernel* mKernel;                   ///< Pointer to kernel function
  //  vectord mTheta;                             ///< Kernel parameters
  const double mRegularizer;  ///< GP likelihood parameters (Normal)

};

/**@}*/

#endif
