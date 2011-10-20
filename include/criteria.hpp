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

#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution


using namespace boost::numeric::ublas;	
using boost::math::factorial;
using boost::math::normal; // typedef provides default type is double.

enum criterium_name{
  expectedImprovement = 1,
  lcb,
  poi,
  gp_hedge,
  greedyAOptimality,
  expectedReturn
};


class Criteria
{
public:

  Criteria()
  {
    n_calls = 0;
    g = 1;     beta = 1;   epsilon = 0;
    g_ei = 0;  g_lcb = 0;  g_poi = 0;
    eta = 100;
  }

  virtual ~Criteria(){}

  inline double evaluate(NonParametricProcess &gp, const vectord &query,
			 criterium_name crit = expectedImprovement)
  {
    n_calls++;

    double yPred, sPred, yMin = gp.getValueAtMinimum(); 
    gp.prediction(query,yPred,sPred);

    switch (crit)
      {
      case expectedImprovement: return negativeExpectedImprovement(yPred,sPred,yMin);
      case lcb: return lowerConfidenceBound(yPred,sPred);
      case poi: return probabilityOfImprovement(yPred,sPred,yMin);
      case gp_hedge: return hedge(yPred,sPred,yMin);
      case greedyAOptimality: return sPred;
      case expectedReturn: return yPred;
      default: std::cout << "Error in criterium" << std::endl; return 0.0;
      }
  }

protected:

  inline double negativeExpectedImprovement(double yPred, double sPred,
					    double yMin)
  {
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
	  sumEI += pow(-1.0,ii)*(fg/(factorial<double>(ii)*factorial<double>(g-ii)))*
	    pow(yNorm,g-ii)*Tact;
	  
	  //roll-up
	  Tm2 = Tm1;   Tm1 = Tact;
	}
      return -1.0 * pow(sPred,g) * sumEI;
    }
         
  }  // negativeExpectedImprovement


  inline double lowerConfidenceBound(double yPred, double sPred)
  {    
    return yPred -  beta*sPred;;
  }

  inline double probabilityOfImprovement(double yPred, double sPred,
					 double yMin)
  {
    return cdf(d,(yMin - yPred + epsilon)/sPred);
  }

  inline double hedge(double yPred, double sPred, double yMin)
  {
    /*    double ei = negativeExpectedImprovement(yPred,sPred,yMin,g);
    double lcb = lowerConfidenceBound(yPred,sPred,beta);
    double poi = probabilityOfImprovement(yPred,sPred,yMin,epsilon);

    double p_ei = exp(eta*g_ei);
    double p_lcb = exp(eta*g_lcb);
    double p_poi = exp(eta*g_poi);
    double sum_p = p_ei + p_lcb + p_poi;

    p_ei /= sum_p; p_lcb /= sum_p; p_poi /= p_poi;
    g_ei += 

    double u = sample();
    if (u < p_ei)
      return ei;
    else if (u < p_lcb+p_ei)
      return lcb;
    else
    return poi;*/
    return 0.0;
  }


  inline void updateCoolingScheme()
  {
    if (n_calls%10)
      g = std::max(1,static_cast<int>(round(g/2.0)));
  }

  
protected:
 
  unsigned int n_calls, g;

  double beta, epsilon, eta; 
  double g_poi,g_lcb,g_ei;

  normal d;


};


#endif
