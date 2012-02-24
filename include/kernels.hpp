/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#ifndef  _KERNELS_HPP_
#define  _KERNELS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include "elementwiseUblas.hpp"
#include "ctypes.h"

namespace kernels
{

  using namespace boost::numeric::ublas;	

  /** 
   * Matern isotropic kernel
   * 
   * @param x1 First point
   * @param x2 Secont point
   * @param param_index 
   *           if < 0 return # of params
   *           if = 0 return kernel value
   *           if > 0 return the derivative of the ith param
   * @param theta kernel length-scale
   * @param order kernel order+1/2
   * 
   * @return see param_index
   */
  inline double MaternIso( const vector<double> &x1, 
			    const vector<double> &x2,
			    size_t param_index,
			    double theta, int order)
  {
    if (param_index < 0)
      return 1; /*Number of hyperparams*/

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
  } 
 
  /** 
   * Gaussian isotropic kernel
   * 
   * @param x1 First point
   * @param x2 Secont point
   * @param param_index 
   *           if < 0 return # of params
   *           if = 0 return kernel value
   *           if > 0 return the derivative of the ith param
   * @param theta kernel length-scale
   * 
   * @return see param_index
   */
  inline double SEIso( const vector<double> &x1, 
		       const vector<double> &x2,
		       size_t param_index,
		       double theta)
  {
    if (param_index < 0)
      return 1; /*Number of hyperparams*/

    double rl = norm_2(x1-x2)/theta;
    double k = rl*rl;
    double result = exp(-k/2);
    double grad = result*k;

    if (param_index == 0) return result;
    else return grad;
  }

 
  /** 
   * Gaussian kernel with Automatic Relevance Determination (ARD) 
   * distance measure
   * 
   * @param x1 First point
   * @param x2 Secont point
   * @param param_index 
   *           if < 0 return # of params
   *           if = 0 return kernel value
   *           if > 0 return the derivative of the ith param
   * @param theta kernel length-scale
   * 
   * @return see param_index
   */
  inline double SEard( const vector<double> &x1, 
		       const vector<double> &x2,
		       size_t param_index,
		       const vector<double> &theta)
  {
    if (param_index < 0)
      return x1.size(); /*Number of hyperparams*/

    vector<double> xd = x1-x2;
    vector<double> ri = ublas_elementwise_div(xd, theta);

    double rl = norm_2(ri);
    double k = rl*rl;
    double result = exp(-k/2);
    double grad = result*sqrt(ri(param_index));

    if (param_index == 0) return result;
    else return grad;
  }


  /** 
   * Common call function for the kernel functions
   * 
   * @param kname name of the kernel function
   * @param x1 First point
   * @param x2 Secont point
   * @param param_index 
   *           if < 0 return # of params
   *           if = 0 return kernel value
   *           if > 0 return the derivative of the ith param
   * @param params for all functions: length-scales 
   *               for Matern: the last component is order+1/2 
   * 
   * @return 
   */
  inline double kernelFunction( kernel_name kname,
				const vector<double> &x1, 
				const vector<double> &x2,
				size_t param_index,
				vector<double> params )
  {
    // TODO: Add asserts for dimension checking depending on kernel function
    switch(kname)
      {
      case k_materniso: return MaternIso(x1,x2,param_index,params(0),params(1));
      case k_seiso:      return SEIso(x1,x2,param_index,params(0));
      case k_seard:      return SEard(x1,x2,param_index,params);
      default: 
	std::cout << "Kernel function not supported!" << std::endl; 
	return 0.0;
      }
  }



}


#endif
