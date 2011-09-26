//
// C++ Implementation: criteria
//
// Description: 
//
//
// Author: Ruben Martinez-Cantin  <rmcantin@unizar.es>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//



#ifndef  _CRITERIA_HPP_
#define  _CRITERIA_HPP_

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;	

namespace criteria
{

  unsigned int factorial(unsigned int no, unsigned int a = 1)
  {
    // termination condition
    if (0 == no || 1 == no)
      return a;
  
    // Tail recursive call
    return factorial(no - 1, no * a);
  } //factorial

  double pdf(double x)
  {
    return (1 / (sqrt(2 * M_PI)) * exp(-(x*x)/2));
  } //pdf
  
  double cdf(double x)
  {
    /** \brief Abromowitz and Stegun approximation of Normal CDF
     * 
     * Extracted from 
     * http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
     * The most used algorithm is algorithm 26.2.17 from Abromowitz and Stegun, 
     * Handbook of Mathematical Functions. It has a maximum absolute error of 7.5e^-8.
     * 
     */
	
    static const double b1 =  0.319381530;
    static const double b2 = -0.356563782;
    static const double b3 =  1.781477937;
    static const double b4 = -1.821255978;
    static const double b5 =  1.330274429;
    static const double p  =  0.2316419;
    static const double c  =  0.39894228;
  
    if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
	      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
    else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
	       ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
  } // cdf



  double negativeExpectedImprovement(double yPred, double sPred,
				     double yMin, int g = 1)
  {
    double yDiff = yMin - yPred; 
    double yNorm = yDiff / sPred;
  
    if (g == 1)
	return -1.0 * ( yDiff * cdf(yNorm) + sPred * pdf(yNorm) );
    else
    {
      double pdfD = pdf(yNorm);
      double Tm2 = cdf(yNorm);
      double Tm1 = pdfD;
      double fg = factorial(g);
      double Tact;
      double sumEI = pow(yNorm,g)*Tm2 - g*pow(yNorm,g-1)*Tm1;

      for (int ii = 2; ii < g; ii++) 
	{
	  Tact = (ii-1)*Tm2 - pdfD*pow(yNorm,ii-1);
	  sumEI += pow(-1.0,ii)*(fg/(factorial(ii)*factorial(g-ii)))*
	    pow(yNorm,g-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      return -1.0 * pow(sPred,g) * sumEI;
    }
         
  }  // negativeExpectedImprovement


  double lowerConfidenceBound(double yPred, double sPred,
			      double beta)
  {    
    return yPred -  beta*sPred;;
  }

  
  double greedyAOptimality(double yPred, double sPred)
  { return sPred; } 

}


#endif
