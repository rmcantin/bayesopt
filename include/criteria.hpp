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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution


using namespace boost::numeric::ublas;	
using boost::math::factorial;
using boost::math::normal; // typedef provides default type is double.

enum criterium_name{
  expectedImprovement = 1,
  lowerConfidenceBound,
  probabilityOfImprovement,
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

  inline double evaluateCriteria(double yPred, double sPred, double yMin, 
				 criterium_name crit = 1)
  {
    switch (crit)
      {
      case expectedImprovement: return negativeExpectedImprovement(yPred,sPred,yMin);
      case lowerConfidenceBound: return lowerConfidenceBound(yPred,sPred);
      case probabilityOfImprovement: return probabilityOfImprovement(yPred,sPred,yMin);
      case greedyAOptimality: sPred;
      case expectedReturn: yPred;
      default: std::cout << "Error in criterium" << std::endl;
      }
  }

  inline double negativeExpectedImprovement(double yPred, double sPred,
					    double yMin, int g = 1)
  {
    double yDiff = yMin - yPred; 
    double yNorm = yDiff / sPred;

    g = 2;
  
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

      for (int ii = 2; ii < g; ii++) 
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


  inline double lowerConfidenceBound(double yPred, double sPred,
				     double beta)
  {    
    return yPred -  beta*sPred;;
  }

  inline double probabilityOfImprovement(double yPred, double sPred,
					 double yMin, double epsilon)
  {
    return cdf(d,(yMin - yPred + epsilon)/sPred);
  }

  inline double hedge(double yPred, double sPred, double yMin, double eta, 
		      double epsilon, double beta, double g = 1)
  {
    double ei = negativeExpectedImprovement(yPred,sPred,yMin,g);
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
      return poi;
  }
  
protected:



  inline int updateCoolingScheme(size_t nTotalIterations,
				 size_t nCurrentIteration)
  {

    double iterPercentaje = static_cast<double>(nTotalIterations) 
      / static_cast<double>(nCurrentIteration);

    if (iterPercentaje < 0.2)
      mG = 20;
    else if (iterPercentaje < 0.3)
      mG = 10; 
    else if (iterPercentaje < 0.4)
      mG = 5;     
    else if (iterPercentaje < 0.5)
      mG = 2;    
    else if (iterPercentaje < 0.7)
      mG = 1;       

    return 1;
  }


 
  const bool mUseCool;
  unsigned int n_calls, g;

  double mLCBparam;                   // LCB = mean - param * std
  criterium_name crit;


  normal d;


}


#endif
