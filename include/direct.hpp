#ifndef _DIRECT_HPP_
#define _DIRECT_HPP_

#include <stdint.h>

#if defined(WIN32) || defined(__MINGW32__) 
  #if defined(DIRect_EXPORTS)
    #define  DllExport __declspec(dllexport)
  #else
    #define  DllExport __declspec(dllimport)
  #endif /* MyLibrary_EXPORTS */
#else /* defined (_WIN32) */
 #define DllExport
#endif


namespace DIRECT
{
	// DIRECT routines as C callable functions
	extern "C" {
	  DllExport int directhead_(
			  int, double *, int *, 
			  double *, double *,
			  double *, int *, int *, int *, int *
			  );

	  int callfcn_   (int *n, double *x, double *f, 
			  int *flag__, int *iidata, 
			  int *iisize, double *ddata, 
			  int *idsize, char *cdata,
			  int *icsize, int cdata_len );	 

	  int negeiwrap_   (int *n, double *x, double *f, 
			  int *flag__, int *iidata, 
			  int *iisize, double *ddata, 
			  int *idsize, char *cdata,
			  int *icsize, int cdata_len );

	  int lcbwrap_   (int *n, double *x, double *f, 
			  int *flag__, int *iidata, 
			  int *iisize, double *ddata, 
			  int *idsize, char *cdata,
			  int *icsize, int cdata_len );

	} // extern "C"

  

	// Type overloads for C++
	inline void direct(int(*fpointer)(int *, double *, double *, 
					  int *, int *,int *, double *,
					  int *, char *, int *, int),
			   double *x, int *n, 
			   double *fmin, double *l,double *u,
			   int *Ierror, int *maxf, int *maxT, void *objPointer)
	{
	  directhead_((uintptr_t)fpointer, x, n, fmin, l, u, Ierror, maxf, maxT, (int*)objPointer);
	}
}

#endif
