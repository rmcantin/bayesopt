#include "direct.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include "krigging.hpp"



 
namespace DIRECT
{

  /* +----------------------------------------------------------+ */
  /* | C wrapper for DIRECT optimization of current criteria    | */
  /* +----------------------------------------------------------+ */
  int evaluate_wrap_   (int *n, double *x, double *f, 
		    int *flag__, int *iidata, 
		    int *iisize, double *ddata, 
		    int *idsize, char *cdata,
		    int *icsize, int cdata_len)
  {
    array_adaptor<double> shared((*n), x);
    vector<double, array_adaptor<double> > sharedN((*n), shared); 
    
    // This is not very clever... but works!
    void *objPointer = iidata;
    InnerOptimization* OPTIMIZER = static_cast<InnerOptimization*>(objPointer);
    
    vector<double> vgrad(n);
    *f = OPTIMIZER->innerEvaluate(sharedN,vgrad);
    *flag__ = 0;
    
    return 0;
  } /* criteriawrap_ */


          /* +------------------------------------------+ */
	  /* | Test function body.                      | */
	  /* | Used to test DIRECT.                     | */
	  /* +------------------------------------------+ */
	  int callfcn_   (int *n, double *x, double *f, 
			  int *flag__, int *iidata, 
			  int *iisize, double *ddata, 
			  int *idsize, char *cdata,
			  int *icsize, int cdata_len)
	  {
	    int i,dim;
	    /* Parameter adjustments */
	    --x;
	    /* Function Body */
	    *f = 10.;
	    dim = *n;
	    for (i = 1; i <= dim; ++i)
	      {
		*f += (x[i] - .53f) * (x[i] - .53f);
	      }
	    *flag__ = 0;
	    return 0;
	  } /* callfcn_ */


}
