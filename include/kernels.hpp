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

namespace kernels
{
  using namespace boost::numeric::ublas;	

  /** 
   * Mattern isotropic kernel
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
  inline double MatternIso( const vector<double> &x1, 
			    const vector<double> &x2,
			    size_t param_index,
			    double theta, int order)
  {
    /** \brief Matern kernel
     * Matern covariance function
     */

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
    double rl = norm_2(x1-x2)/theta;
    double k = rl*rl;
    double result = exp(-k/2);
    double grad = result*k;

    if (param_index == 0) return result;
    else return grad;
  }

}


#endif
