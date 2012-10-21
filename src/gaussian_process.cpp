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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/normal.hpp> 
#include "gaussian_process.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"

using boost::numeric::ublas::inplace_solve;
using boost::numeric::ublas::lower_tag;
using boost::numeric::ublas::lower;
using boost::numeric::ublas::inner_prod;

GaussianProcess::GaussianProcess( double noise ):
  NonParametricProcess(noise)
{}  // Constructor


GaussianProcess::~GaussianProcess()
{} // Default destructor



double GaussianProcess::negativeLogLikelihood(size_t index)
{
  matrixd K = computeCorrMatrix();
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  // Compute the likelihood
  vectord alpha(mGPY);
  cholesky_solve(L,alpha,lower());
  double loglik = .5*inner_prod(mGPY,alpha) + trace(L) + n*0.91893853320467; 
  // 0.9183... = log(2*pi)/2

  return loglik;
}


int GaussianProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  double kn = (*mKernel)(query, query);
  vectord kStar = computeCrossCorrelation(query);

  return predictionChol(kn,kStar,yPred,sPred);
}

int GaussianProcess::predictionChol(double kn, const vectord& kStar,
				    double& yPred, double& sPred)
{
  yPred = inner_prod(kStar,mAlphaV);

  vectord vd(kStar);
  inplace_solve(mL,vd,lower_tag());
  sPred = sqrt(kn - inner_prod(vd,vd));
  
  return 1;
}

int GaussianProcess::predictionInv(double kn, const vectord& kStar,
				    double& yPred, double& sPred)
{
  vectord rInvR = prod(kStar,mInvR);
  yPred = inner_prod(rInvR,mGPY);
  sPred = sqrt(kn - inner_prod(rInvR,kStar));

  return 1;
}

int GaussianProcess::precomputePrediction()
{ return precomputeChol(); }

int GaussianProcess::precomputeChol()
{
  size_t n = mGPY.size();
  
  mAlphaV.resize(n,false);
  mAlphaV = mGPY;
  cholesky_solve(mL,mAlphaV,lower());

  return 1; 
}

int GaussianProcess::precomputeInv()
{ return 1;}
	
double GaussianProcess::negativeExpectedImprovement(const vectord& query,
						    size_t g)
{
  double yPred,sPred;
  double yMin = getValueAtMinimum();
  prediction(query,yPred,sPred);
  
  using boost::math::factorial;

  boost::math::normal d;
  double yDiff = yMin - yPred; 
  double yNorm = yDiff / sPred;
  
  if (g == 1)
    return -1.0 * ( yDiff * cdf(d,yNorm) + sPred * pdf(d,yNorm) );
  else
    {
      double pdfD = pdf(d,yNorm);
      double Tm2 = cdf(d,yNorm);
      double Tm1 = pdfD;
      double fg = factorial<double>(g);
      double Tact;
      double sumEI = pow(yNorm,g)*Tm2 - g*pow(yNorm,g-1)*Tm1;

      for (unsigned int ii = 2; ii < g; ii++) 
	{
	  Tact = (ii-1)*Tm2 - pdfD*pow(yNorm,ii-1);
	  sumEI += pow(-1.0,ii)* 
	    (fg / ( factorial<double>(ii)*factorial<double>(g-ii) ) )*
	    pow(yNorm,g-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      return -1.0 * pow(sPred,g) * sumEI;
    }
  
}  // negativeExpectedImprovement

double GaussianProcess::lowerConfidenceBound(const vectord& query,
					     double beta)
{    
  double yPred,sPred;
  prediction(query,yPred,sPred);
  return yPred - beta*sPred;;
}  // lowerConfidenceBound


double GaussianProcess::negativeProbabilityOfImprovement(const vectord& query,
							 double epsilon)
{
  boost::math::normal d;
  double yMin = getValueAtMinimum();
  double yPred,sPred;
  prediction(query,yPred,sPred);
  return -cdf(d,(yMin - yPred + epsilon)/sPred);
}  // negativeProbabilityOfImprovement
