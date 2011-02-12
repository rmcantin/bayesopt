#ifndef _NLOPTWPR_HPP_
#define _NLOPTWPR_HPP_


namespace NLOPT_WPR
{
  extern "C" {

  /** 
   * \brief Negative Expected Improvement C wrapper for NLOPT
   * 
   * @param n size of query
   * @param x query point
   * @param grad gradient if available (not used now)
   * @param my_func_data any other function data
   * 
   * @return negative EI 
   */
  double negeiwrap_nlopt (unsigned int n, const double *x,
			  double *grad, void *my_func_data);
  
  /** 
   * \brief Lower Confidence Bound C wrapper for NLOPT
   * 
   * @param n size of query
   * @param x query point
   * @param grad gradient if available (not used now)
   * @param my_func_data any other function data
   * 
   * @return LCB
   */
  double lcbwrap_nlopt (unsigned int n, const double *x,
			double *grad, void *my_func_data);

  }
}

#endif
