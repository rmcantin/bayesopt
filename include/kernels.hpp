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
  inline double MatternIso( const vector<double> &x1, 
			    const vector<double> &x2,
			    size_t param_index,
			    double theta, int order)
  {
    /** \brief Matern kernel
     * Matern covariance function
     */
    double k,grad;

    double r = sqrt(order) * norm_2(x1-x2)/theta;
    double er = exp(-r);

    switch (order)
      {
      case 1: grad = r*er;  k=er; break;
      case 3: grad = r*r*er; k=(1+r)*er; break;
      case 5: grad = r*(1+r)/3*r*er; k=(1+r*(1+r/3))*er; break;
      default: 
	std::cout << "Error: not suported kernel." << std::endl;
	return 0.0;
      }
    if (param_index == 0) return k;
    else return grad;
  }  // correlationFunction
  
  inline double SEIso( const vector<double> &x1, 
		       const vector<double> &x2,
		       size_t param_index,
		       double theta)
  {
    double rl = norm_2(x1-x2)/theta;
    double k = rl*rl;
    double result = exp(-k/2);
    double grad = result*k;

    if (param_index == 0) return result;
    else return grad;
  }

}


#endif
