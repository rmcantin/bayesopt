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
#ifndef _DIRECT_HPP_
#define _DIRECT_HPP_

#include <stdint.h>

// DLL magic for windows compilers
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

	  /** 
	   * Wrapper of inner optimization to be evaluated by DIRECT
	   * 
	   * @param n # of dimensions
	   * @param x input point
	   * @param grad returns gradient evaluation
	   * @param my_func_data pointer to the InnerOptimization object
	   * 
	   * @return function evaluation
	   */  
	  int evaluate_wrap_   (int *n, double *x, double *f, 
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
