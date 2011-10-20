

#include "basicgaussprocess.hpp"
#include "cholesky.hpp"
#include "trace.hpp"
  

BasicGaussianProcess::BasicGaussianProcess():
  NonParametricProcess(), mTheta(KERNEL_THETA),
  mRegularizer(DEF_REGULARIZER)
{} // Default constructor

BasicGaussianProcess::BasicGaussianProcess( double theta, 
					    double noise):
  NonParametricProcess(), mTheta(theta),
  mRegularizer(noise)
{}  // Constructor


BasicGaussianProcess::~BasicGaussianProcess()
{} // Default destructor



double BasicGaussianProcess::negativeLogLikelihood(double param,
						   double &grad)
{
  matrixd K = computeCorrMatrix(mRegularizer,0);
  matrixd dK = computeCorrMatrix(0,1);
  size_t n = K.size1();
  matrixd L(n,n);
  cholesky_decompose(K,L);
  
  vectord alpha(mGPY);
  matrixd inverse = eyed(n);

  boost::numeric::ublas::inplace_solve(L,alpha,lower_tag());
  boost::numeric::ublas::inplace_solve(L,inverse,lower_tag());

  double loglik = -.5*inner_prod(mGPY,alpha) - trace(L) - n*0.91893853320467; //log(2*pi)/2
  grad = 0.5 * trace_prod(outer_prod(alpha,alpha) - inverse, dK);

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
	

int BasicGaussianProcess::fitGP()
{
  size_t nSamples = mGPXX.size();
  for (size_t ii=0; ii<nSamples; ii++)
    checkBoundsY(ii);

  int error = computeInverseCorrMatrix(mRegularizer);

  if (error < 0)
    return error;

  return 1;
} // fitGP


int BasicGaussianProcess::addNewPointToGP(const vectord &Xnew, 
					  double Ynew)
{
  size_t nSamples = mGPXX.size();
  size_t XDim = mGPXX[1].size();
  size_t NewDim = Xnew.size();
  
  svectord colU(nSamples+1,1.0);
  vectord Li(nSamples);
  vectord wInvR(nSamples);
  double wInvRw;
  double selfCorrelation, Ni;
  
  if (XDim != NewDim)
    {
      std::cout << "Dimensional Error" << std::endl;
      return -1;
    }
    
  vectord correlationNewValue = computeCrossCorrelation(Xnew);
  
  selfCorrelation = correlationFunction(Xnew, Xnew) + mRegularizer;
  
  noalias(wInvR) = prod(correlationNewValue,mInvR);
  wInvRw = inner_prod(wInvR,correlationNewValue);
  Ni = 1/(selfCorrelation + wInvRw);
  noalias(Li) = -Ni * wInvR;
  mInvR += outer_prod(Li,Li) / Ni;
  
  //TODO: There must be a better way to do this.
  mInvR.resize(nSamples+1,nSamples+1);
  
  Li.resize(nSamples+1);
  Li(nSamples) = Ni;
  
  row(mInvR,nSamples) = Li;
  column(mInvR,nSamples) = Li;

  addSample(Xnew,Ynew);
  checkBoundsY(nSamples);
  
  return 1;
} // addNewPointToGP

