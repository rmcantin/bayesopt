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

#ifndef  _BASICGAUSSPROCESS_HPP_
#define  _BASICGAUSSPROCESS_HPP_

#include "krigwpr.h"
#include "nonparametricprocess.hpp"
#include "kernels.hpp"
#include "meanfuncs.hpp"
	
 
/** \addtogroup BayesOptimization */
/*@{*/


class BasicGaussianProcess: public NonParametricProcess
{
public:
  BasicGaussianProcess();  
  BasicGaussianProcess( double theta, double noise );
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
  inline double correlationFunction( const vectord &x1,const vectord &x2,
				     size_t param_index = 0 )
  { return kernels::MatternIso(x1,x2,param_index,mTheta,3); }

  inline double meanFunction( const vectord &x )
  { return means::Zero(x); }


  int precomputeGPParams();


protected:
  double mTheta;                      // Kernel parameters
  const double mDelta2,mRegularizer;  // GP prior parameters (Normal)

  // Precomputed GP prediction operations
  vectord mUInvR;              
  double mUInvRUDelta;
  vectord mYUmu;

};


/**@}*/
// end namespaces

#endif
