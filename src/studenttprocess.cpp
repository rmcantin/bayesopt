#include <boost/math/distributions/students_t.hpp> // for student t distribution

#include "studenttprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"

  
StudentTProcess::StudentTProcess(double noise):
  NonParametricProcess(noise)
{}  // Constructor



StudentTProcess::~StudentTProcess()
{} // Default destructor




double StudentTProcess::negativeLogLikelihood(size_t index)
{
  matrixd K = computeCorrMatrix(0);
  size_t n = K.size1();
  
  matrixd L(n,n);
  cholesky_decompose(K,L);

  vectord colU(n);

  //TODO: Replace by transform
  for (size_t ii=0; ii< n; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

  vectord alphU(colU);
  boost::numeric::ublas::inplace_solve(L,alphU,boost::numeric::ublas::lower_tag());
  double eta = inner_prod(colU,alphU);
  
  vectord alphY(mGPY);
  boost::numeric::ublas::inplace_solve(L,alphY,boost::numeric::ublas::lower_tag());
  double mu     = inner_prod(colU,alphY) / eta;
  double YInvRY = inner_prod(mGPY,alphY);
    
  double sigma = (YInvRY - mu*mu*eta) / (n-1);

  double negloglik = 0.5*( (n-1)*log(sigma) + trace(L) + log(eta) );

  return negloglik;
}


int StudentTProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  size_t n = mGPXX.size();
  vectord rInvR(n);
  double kn;
  double uInvRr, rInvRr;
  double meanf = meanFunction(query);
  
  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  svectord colMu(n,mMu);
  vectord yumu = mGPY - meanf*colMu;
  
  yPred = meanf*mMu + inner_prod( rInvR, yumu );
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			/ mUInvRUDelta ) );

  return n-1;
}
	

int StudentTProcess::precomputeGPParams()
{
  size_t nSamples = mGPXX.size();
  vectord colU(nSamples);

  //TODO: Replace by transform
  for (size_t ii=0; ii< nSamples; ii++) 
    colU(ii) = meanFunction(mGPXX[ii]);

  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU);
  
  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  
  mSig = (YInvRY - mMu*mMu*mUInvRUDelta) / (nSamples-1);

  return 1;
}


	
double StudentTProcess::negativeExpectedImprovement(double yPred, double sPred,
						    double yMin, size_t g)
{
  size_t dof = mGPXX.size() - 1;
  boost::math::students_t d(dof);

  double yDiff = yMin - yPred; 
  double yNorm = yDiff / sPred;
  
  if (g != 1)
    {
      std::cout << "Students t EI with exponent not yet supported." << std::endl;
      return 0.0;
    }
  else
    {
      return -1.0 * ( yDiff * cdf(d,yNorm) + (dof*sPred+yNorm*yDiff)/(dof-1) * pdf(d,yNorm) );
    }
  
}  // negativeExpectedImprovement

double StudentTProcess::lowerConfidenceBound(double yPred, double sPred,
					     double beta)
{    
  size_t n = mGPXX.size();
  return yPred - beta*sPred/sqrt(n);
}  // lowerConfidenceBound

double StudentTProcess::negativeProbabilityOfImprovement(double yPred, double sPred,
							      double yMin, double epsilon)
{  
  size_t dof = mGPXX.size() - 1;
  boost::math::students_t d(dof);
  return -cdf(d,(yMin - yPred + epsilon)/sPred);
}  // negativeProbabilityOfImprovement
