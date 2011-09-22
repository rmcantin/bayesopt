//
// C++ Implementation: kernels
//
// Description: 
//
//
// Author: Ruben Martinez-Cantin  <rmcantin@unizar.es>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//



#ifndef  _KERNELS_HPP_
#define  _KERNELS_HPP_

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;	

namespace kernels
{
  double MatternIso( const vector<double> &x1, 
		     const vector<double> &x2,
		     double &grad,
		     double theta, int order)
  {
    /** \brief Matern kernel
     * Matern covariance function
     */

    double r = sqrt(order) * norm_2(x1-x2)/theta;
    double er = exp(-r);

    switch (order)
      {
      case 1: grad = r*er;  return er;
      case 3: grad = r*r*er; return (1+r)*er;
      case 5: grad = r*(1+r)/3*r*er; return (1+r*(1+r/3))*er;
      defaul: 
	std::cout << "Error: not suported kernel." << std::endl;
	return 1;
      }

  }  // correlationFunction
  
  double SEIso( const vector<double> &x1, 
		const vector<double> &x2,
		double &grad,
		double theta)
  {
    double rl = norm_2(x1-x2)/theta;
    double k = rl*rl;
    double result = exp(-k/2);
    
    grad = result*k;

    return result;
  }

}


#endif
