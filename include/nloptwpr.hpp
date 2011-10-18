#ifndef _NLOPTWPR_HPP_
#define _NLOPTWPR_HPP_


namespace NLOPT_WPR
{
  extern "C" {

  /** 
   * \brief C wrapper for NLOPT to evaluate the criteria for a query
   * 
   * @param n size of query
   * @param x query point
   * @param grad gradient if available (not used now)
   * @param my_func_data any other function data
   * 
   * @return negative EI 
   */
  double evaluate_nlopt (unsigned int n, const double *x,
				  double *grad, void *my_func_data);
  
  }
}

#endif
