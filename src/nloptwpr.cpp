#include "nloptwpr.hpp"
// BOOST Libraries
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "krigging.hpp"

using namespace boost::numeric::ublas;	

namespace NLOPT_WPR
{
/* +-------------------------------------------------------+ */
/* | Negative Expected Improvement C wrapper for NLOPT     | */
/* +-------------------------------------------------------+ */
  double evaluate_nlopt (unsigned int n, const double *x,
				  double *grad, void *my_func_data)

  {
    double xcopy[128];
    for (unsigned int i=0;i<n;i++)
      xcopy[i] = x[i];
    array_adaptor<double> shared(n, xcopy);
    vector<double, array_adaptor<double> > sharedN(n, shared); 
    
    // This is not very clever... but works!
    void *objPointer = my_func_data;
    InnerOptimization* GAUSSIAN_PROCESS = static_cast<InnerOptimization*>(objPointer);
    
    double f =  GAUSSIAN_PROCESS->innerEvaluate(sharedN);
    
    return f;
  } /* evaluate_criteria_nlopt */

}
