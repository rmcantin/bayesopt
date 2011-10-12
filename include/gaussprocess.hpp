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

#include "krigwpr.h"
#include "nonparametricprocess.hpp"
	
 
/** \addtogroup BayesOptimization */
/*@{*/


class GaussianProcess: public NonParametricProcess
{
public:
  GaussianProcess();  
  GaussianProcess( double theta, double p,
		   double alpha, double beta, 
		   double delta, double noise);
  GaussianProcess( gp_params params );
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

  //  int marginalLikelihood(const vectord &query,
			 
			 

  /** Computes the GP based on mGPXX
   *  This function is hightly inefficient O(N^3). Use it only at 
   *  the begining
   */
  int fitGP();



  /** Add new point efficiently using Matrix Decomposition Lemma
   *  for the inversion of the correlation matrix. Maybe it is faster
   *  to just construct and invert a new matrix each time.
   */   
  int addNewPointToGP( const vectord &Xnew,
		       double Ynew);

  inline double getTheta()
  { return mTheta; }

  inline void setTheta( double theta )
  { mTheta = theta; }


protected:
  double correlationFunction( const vectord &x1,
			      const vectord &x2 );


  double correlationFunction( const vectord &x1, 
			      const vectord &x2,
			      double param, double &grad);

  int precomputeGPParams();

  double meanFunction( const vectord &x);

protected:
  double mTheta, mP;                  // Kernel parameters
  const double mAlpha, mBeta;         // GP prior parameters (Inv-Gamma)
  const double mDelta2,mRegularizer;  // GP prior parameters (Normal)

  double mMu, mSig;                   // GP posterior parameters

  // Precomputed GP prediction operations
  vectord mUInvR;              
  double mUInvRUDelta;
  vectord mYUmu;

};


/**@}*/
// end namespaces

#endif
