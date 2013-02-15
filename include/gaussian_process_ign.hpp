/** \file gaussian_process_ign.hpp 
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

#ifndef  _GAUSSIAN_PROCESS_IGN_HPP_
#define  _GAUSSIAN_PROCESS_IGN_HPP_

#include "nonparametricprocess.hpp"
 
/** \addtogroup  NonParametricProcesses */
/**@{*/

/**
 * \brief Gaussian process with normal-inverse-gamma hyperprior 
 *        on mean and signal variance parameters.
 */
class GaussianProcessIGN: public NonParametricProcess 
{
public:
  GaussianProcessIGN(size_t dim, double noise, double alpha,
		     double beta, double delta);

  virtual ~GaussianProcessIGN();

  /** 
   * \brief Function that returns the prediction of the GP for a query point
   * in the hypercube [0,1].
   * 
   * @param query in the hypercube [0,1] to evaluate the Gaussian process
   * @return pointer to the probability distribution.
   */	
  ProbabilityDistribution* prediction(const vectord &query);

protected:

  /** 
   * \brief Computes the negative log likelihood and its gradient of the data.
   * 
   *
   * @return value negative log likelihood
   */
  double negativeLogLikelihood();

  /** 
   * \brief Precompute some values of the prediction that do not depends on
   * the query
   * @return error code
   */
  int precomputePrediction();

private:
  const double mAlpha, mBeta;         //!< GP prior parameters (Inv-Gamma)
  const double mDelta2;               //!< GP prior parameters (Normal)

  double mMu, mSig;                   //!< GP posterior parameters

  //! Precomputed GP prediction operations
  vectord mUInvR;
  vectord mInvRy;
  double mUInvRUDelta;
  vectord mAlphaV;

};
/**@}*/

#endif
