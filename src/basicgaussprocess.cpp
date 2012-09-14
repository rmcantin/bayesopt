#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution

#include "basicgaussprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"


BasicGaussianProcess::BasicGaussianProcess( double noise ):
  NonParametricProcess(noise)
{}  // Constructor


BasicGaussianProcess::~BasicGaussianProcess()
{} // Default destructor



double BasicGaussianProcess::negativeLogLikelihood(size_t index)
{
  matrixd K = computeCorrMatrix();
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);

  // Compute the likelihood
  vectord alpha(mGPY);
  boost::numeric::ublas::inplace_solve(L,alpha,boost::numeric::ublas::lower_tag());
  double loglik = .5*inner_prod(mGPY,alpha) + trace(L) + n*0.91893853320467; //log(2*pi)/2

  return loglik;
}


int BasicGaussianProcess::prediction( const vectord &query,
				      double& yPred, double& sPred)
{
  vectord rInvR(mGPXX.size());
  double kn;
  double rInvRr;

  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
    
  yPred = inner_prod( rInvR, mGPY );
  sPred = sqrt(kn - rInvRr);

  return 1;
}
	
double BasicGaussianProcess::negativeExpectedImprovement(double yPred, 
							 double sPred,
							 double yMin, 
							 size_t g)
{
  
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

double BasicGaussianProcess::lowerConfidenceBound(double yPred, double sPred,
						  double beta)
{    
  return yPred - beta*sPred;;
}  // lowerConfidenceBound

double BasicGaussianProcess::negativeProbabilityOfImprovement(double yPred, 
							      double sPred,
							      double yMin, 
							      double epsilon)
{
  boost::math::normal d;
  return -cdf(d,(yMin - yPred + epsilon)/sPred);
}  // negativeProbabilityOfImprovement
