/** -*- c++ -*- \file bayesoptbase.hpp \brief Bayesian optimization module */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/


#ifndef  _DLL_STUFF_HPP_
#define  _DLL_STUFF_HPP_

/* use stdcall convention under Windows, since this seems to
   be more standard there and is important for calling from .NET */
#if defined(_WIN32) || defined(__WIN32__)
#  if defined(__GNUC__)
#    define BAYESOPT_STDCALL __attribute__((stdcall))
#  elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#    define BAYESOPT_STDCALL __stdcall
#  else
#    define BAYESOPT_STDCALL
#  endif
#else
#  define BAYESOPT_STDCALL
#endif

/* WINDOWS DLLs stuff */
#if defined (BAYESOPT_DLL) && (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
  #if defined(bayesopt_EXPORTS)
    #define  BAYESOPT_API __declspec(dllexport)
	#define  BAYESOPT_C_API(T) extern __declspec(dllexport) T BAYESOPT_STDCALL
  #else
    #define  BAYESOPT_API __declspec(dllimport)
	#define  BAYESOPT_C_API(T) extern __declspec(dllimport) T BAYESOPT_STDCALL
  #endif /* MyLibrary_EXPORTS */
#else /* defined (_WIN32) */
 #define BAYESOPT_API
 #define BAYESOPT_C_API(T) extern T BAYESOPT_STDCALL
#endif

        
        
#endif
