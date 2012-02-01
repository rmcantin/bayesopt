/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#ifndef  _STUDENT_T_PROCESS_HPP_
#define  _STUDENT_T_PROCESS_HPP_

#include "nonparametricprocess.hpp"
#include "kernels.hpp"
#include "meanfuncs.hpp"
 
/** \addtogroup BayesOptimization */
/*@{*/


/**
 * \brief Student's t process with hyperpriors on mean and signal variance parameters.
 */
class StudentTProcess: public NonParametricProcess
{
public:
  StudentTProcess( double theta = KERNEL_THETA, 
		   double noise = DEF_REGULARIZER);
  virtual ~StudentTProcess();

  /** 
   * Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query point in the hypercube [0,1] to evaluate the Gaussian process
   * @param yPred mean of the predicted Gaussian distribution
   * @param sPred std of the predicted Gaussian distribution
   * 
   * @return if positive: degrees of freedom, if negative: error code.
   */	
  int prediction(const vectord &query,
  		 double& yPred, double& sPred);


  /** 
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  double negativeLogLikelihood(double& grad,
			       size_t index = 1);			 

  
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
  double negativeExpectedImprovement(double yPred, double sPred,
				     double yMin, size_t g = 1);


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
  double lowerConfidenceBound(double yPred, double sPred,
			      double beta = 1);

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
  double negativeProbabilityOfImprovement(double yPred, double sPred,
					  double yMin, double epsilon = 0.1);

  /** Function to sample from the generalized Student's t distribution that appears
   *  in Pattern Recognition and Machine Learning from C.M. Bishop
   */
  double sample_query(const vectord& query, 
		      randEngine& eng)
  { 
    double y,s;
    size_t n = mGPXX.size() - 1;
    prediction(query,y,s);
    randNFloat normal(eng,normalDist(y,s));
    randGFloat gamma(eng,gammaDist(n/2));
    return normal() / sqrt(2*gamma()/n);
  }


protected:
  inline double correlationFunction( const vectord &x1, const vectord &x2,
				     size_t param_index = 0 )
  { 
    vectord th = mTheta;    th.resize(th.size()+1);   th(th.size()-1) = 3;
    return kernels::kernelFunction(k_matterniso,x1,x2,param_index,th); 
  }

  inline double meanFunction( const vectord &x )
  { return means::One(x); }


  int precomputeGPParams();


protected:

  double mMu, mSig;                   // GP posterior parameters

  // Precomputed GP prediction operations
  vectord mUInvR;              
  double mUInvRUDelta;
};


/**@}*/
// end namespaces

#endif
