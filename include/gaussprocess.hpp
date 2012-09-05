/** \file gaussprocess.hpp 
    \brief Gaussian process with normal-inverse-gamma hyperprior 
           on mean and signal variance parameters. */
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

#ifndef  _GAUSSPROCESS_HPP_
#define  _GAUSSPROCESS_HPP_

//#include "nonparametricprocess.hpp"
#include "basicgaussprocess.hpp"
 
/** \addtogroup  NonParametricProcesses */
/**@{*/


/**
 * \brief Gaussian process with normal-inverse-gamma hyperprior 
 *        on mean and signal variance parameters.
 */
class GaussianProcess: public BasicGaussianProcess 
{
public:
  GaussianProcess( double noise = DEF_REGULARIZER,
		   double alpha = PRIOR_ALPHA, 
		   double beta  = PRIOR_BETA, 
		   double delta = PRIOR_DELTA_SQ );

  virtual ~GaussianProcess();

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
	
  /** 
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  double negativeLogLikelihood(size_t index = 1);			 
			 		 

protected:

  int precomputeGPParams();


protected:
  const double mAlpha, mBeta;         // GP prior parameters (Inv-Gamma)
  const double mDelta2;               // GP prior parameters (Normal)

  double mMu, mSig;                   // GP posterior parameters

  // Precomputed GP prediction operations
  vectord mUInvR;              
  double mUInvRUDelta;

};


/**@}*/
// end namespaces

#endif
