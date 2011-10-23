#include "basicgaussprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"
  

BasicGaussianProcess::BasicGaussianProcess( double theta, 
					    double noise):
  NonParametricProcess(theta,noise)
{
  setAlgorithm(bobyqa);
  setLimits(0.,100.);
}  // Constructor


BasicGaussianProcess::~BasicGaussianProcess()
{} // Default destructor



double BasicGaussianProcess::negativeLogLikelihood(double &grad,
						   size_t index)
{
  matrixd K = computeCorrMatrix(mRegularizer,0);
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  // Compute the likelihood
  vectord alpha(mGPY);
  boost::numeric::ublas::inplace_solve(L,alpha,boost::numeric::ublas::lower_tag());
  double loglik = .5*inner_prod(mGPY,alpha) + trace(L) + n*0.91893853320467; //log(2*pi)/2

#if 0 // BUG:
  // Compute the ith derivative
  if (index > 0)
    {
      matrixd inverse = eyed(n);
      boost::numeric::ublas::inplace_solve(L,inverse,lower_tag());
      matrixd dK = computeCorrMatrix(0,index);
      grad = -.5 * trace_prod(outer_prod(alpha,alpha) - inverse, dK);
    }
#endif

  return loglik;
}


int BasicGaussianProcess::prediction( const vectord &query,
				      double& yPred, double& sPred)
{
  vectord rInvR(mGPXX.size());
  double kn;
  double rInvRr;

  vectord colR = computeCrossCorrelation(query);
  kn = correlationFunction(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
    
  yPred = inner_prod( rInvR, mGPY );
  sPred = sqrt(kn - rInvRr);

  return 1;
}
	
