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

#include "specialtypes.hpp"
#include "cholesky.hpp"
#include "inneroptimization.hpp"	


/**
 * \brief Abstract class to implement non-parametric processes
 */
class NonParametricProcess: public InnerOptimization
{
public:
  NonParametricProcess(double theta = KERNEL_THETA,
		       double noise = DEF_REGULARIZER):
    InnerOptimization(),  mTheta(theta), mRegularizer(noise)
  { mMinIndex = 0; mMaxIndex = 0; }
  
  virtual ~NonParametricProcess(){ }

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
			 double& yPred, double& sPred)
  {return 1;}  

  /** 
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  virtual double negativeLogLikelihood(double& grad,
				       size_t index = 1)
  {return 0.0;}
 
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
  {return 0.0;}

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
  {return 0.0;}

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
						  double yMin, double epsilon = 0.1)
  {return 0.0;}


		 		 
  /** 
   *  Computes the GP based on mGPXX
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the begining
   * 
   * @return error code
   */
  int fitGP()
  {
    size_t nSamples = mGPXX.size();
    for (size_t ii=0; ii<nSamples; ii++)
      checkBoundsY(ii);
  
    vectord th = svectord(1,mTheta);  

    std::cout << "Initial theta: " << mTheta << " "<<th.size()<< std::endl;
    innerOptimize(th);
    setTheta(th(0));
    std::cout << "Final theta: " << mTheta << std::endl;

    int error = computeInverseCorrMatrix(mRegularizer);

    if (error < 0)
      return error;

    return precomputeGPParams();
  } // fitGP

  
  /** 
   *  Add new point efficiently using Matrix Decomposition Lemma
   *  for the inversion of the correlation matrix. Maybe it is faster
   *  to just construct and invert a new matrix each time.
   * 
   * @return error code
   */   
  int addNewPointToGP( const vectord &Xnew,
		       double Ynew)
  {
    size_t nSamples = mGPXX.size();
    size_t XDim = mGPXX[1].size();
  
    vectord Li(nSamples);
    vectord wInvR(nSamples);
    double wInvRw;
    double selfCorrelation, Ni;
  
    if (XDim != Xnew.size())
      {
	std::cout << "Dimensional Error" << std::endl;
	return -1;
      }
    
    vectord correlationNewValue = computeCrossCorrelation(Xnew);
  
    selfCorrelation = correlationFunction(Xnew, Xnew) + mRegularizer;
  
    noalias(wInvR) = prod(correlationNewValue,mInvR);
    wInvRw = inner_prod(wInvR,correlationNewValue);
    Ni = 1/(selfCorrelation + wInvRw);
    noalias(Li) = -Ni * wInvR;
    mInvR += outer_prod(Li,Li) / Ni;
  
    //TODO: There must be a better way to do this.
    mInvR.resize(nSamples+1,nSamples+1);
  
    Li.resize(nSamples+1);
    Li(nSamples) = Ni;
  
    row(mInvR,nSamples) = Li;
    column(mInvR,nSamples) = Li;

    addSample(Xnew,Ynew);
    checkBoundsY(nSamples);
  
    return precomputeGPParams();
  } // addNewPointToGP


  // Getters and setters
  inline void setSamples(matrixd x, vectord y)
  {
    size_t nPoints = x.size1();
    
    for (size_t i=0; i<nPoints; ++i)
      mGPXX.push_back(row(x,i));

    mGPY = y;
  }

  inline void addSample(vectord x, double y)
  {
    mGPXX.push_back(x);
    mGPY.resize(mGPY.size()+1);  mGPY(mGPY.size()-1) = y;
  }

  inline vectord getPointAtMinimum()
  { return mGPXX[mMinIndex]; }

  inline double getValueAtMinimum()
  { return mGPY(mMinIndex); }

  inline double getTheta()
  { return mTheta; }

  inline void setTheta( double theta )
  { mTheta = theta; }

  inline double innerEvaluate(const vectord& query, 
			      vectord& grad)
  { 
    setTheta(query(0));
    return negativeLogLikelihood(grad(0),1);
  }

protected:
  virtual double correlationFunction( const vectord &x1,const vectord &x2, 
				      size_t param_index = 0)
  {return 0.0;}

  virtual double meanFunction( const vectord &x, size_t param_index = 0 )  
  {return 0.0;}

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
  int computeInverseCorrMatrix( double noise )
  {
    size_t nSamples = mGPXX.size();
    if ( (nSamples != mInvR.size1()) || (nSamples != mInvR.size2()) )
      mInvR.resize(nSamples,nSamples);
    
    matrixd corrMatrix = computeCorrMatrix(noise);

    //return InvertMatrix(corrMatrix,mInvR);
    return inverse_cholesky(corrMatrix,mInvR);
  }

  matrixd computeCorrMatrix( double noise , size_t dth_index = 0)
  {
    size_t nSamples = mGPXX.size();
    matrixd corrMatrix(nSamples,nSamples);
  
    for (size_t ii=0; ii< nSamples; ii++)
      {
	for (size_t jj=0; jj < ii; jj++)
	  {
	    corrMatrix(ii,jj) = correlationFunction(mGPXX[ii], mGPXX[jj], dth_index);
	    corrMatrix(jj,ii) = corrMatrix(ii,jj);
	  }
	corrMatrix(ii,ii) = correlationFunction(mGPXX[ii],mGPXX[ii], dth_index);
	if (dth_index == 0) 
	  corrMatrix(ii,ii) += noise;
      }
    return corrMatrix;
  }


  inline vectord computeCrossCorrelation(const vectord &query)
  {
    vectord knx(mGPXX.size());

    //TODO: Replace by transform
    for (size_t ii=0; ii<mGPXX.size(); ++ii)
	knx(ii) = correlationFunction(mGPXX[ii], query);

    return knx;
  }


  inline void checkBoundsY( size_t i )
  {
    if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
    else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
  };

protected:
  vecOfvec mGPXX;                     // TODO:Data inputs
  vectord mGPY;                       // Data values
	
  // Precomputed GP prediction operations
  covMatrix mInvR;                   // Inverse Correlation matrix

  size_t mMinIndex, mMaxIndex;

  double mTheta;                      // Kernel parameters
  const double mRegularizer;  // GP prior parameters (Normal)

};

#endif
