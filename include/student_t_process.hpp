/** \file studenttprocess.hpp
    \brief Student's t process with Jeffreys hyperprior 
           on mean and signal variance parameters. */
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


#ifndef  _STUDENT_T_PROCESS_HPP_
#define  _STUDENT_T_PROCESS_HPP_

#include "nonparametricprocess.hpp"
 
/** \addtogroup  NonParametricProcesses */
/**@{*/


/**
 * \brief Student's t process with Jeffreys hyperprior 
 *        on mean and signal variance parameters.
 */
class StudentTProcess: public NonParametricProcess
{
public:
  StudentTProcess(double noise = DEFAULT_NOISE);
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
   * \brief Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query in the hypercube [0,1] to evaluate the Gaussian process
   * @return pointer to the probability distribution.
   */	
  ProbabilityDistribution* prediction(const vectord &query);


  /** 
   * Computes the negative log likelihood and its gradient of the data.
   * 
   * @param grad gradient of the negative Log Likelihood
   * @param param value of the param to be optimized
   * 
   * @return value negative log likelihood
   */
  double negativeLogLikelihood();			 

  


protected:

  int precomputePrediction();


protected:

  double mMu, mSig;                   // GP posterior parameters

  // Precomputed GP prediction operations
  vectord mUInvR;     
  vectord mInvRy;         
  double mUInvRUDelta;
};


/**@}*/


#endif
