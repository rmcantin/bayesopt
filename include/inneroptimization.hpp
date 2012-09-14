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

#ifndef __INNEROPTIMIZATION_HPP__
#define __INNEROPTIMIZATION_HPP__

#include "specialtypes.hpp"

// We plan to add more in the future since nlopt actually support many of them
enum innerOptAlgorithms {
  DIRECT, // Global optimization
  LBFGS, // Local, derivative based
  BOBYQA, // Local, derivative free
  COMBINED // Global exploration, local refinement
};


class InnerOptimization
{
public:
  InnerOptimization()
  { 
    alg = DIRECT;    mDown = 0.;    mUp = 1.;
  };

  virtual ~InnerOptimization(){};

  /** 
   * Set the optimization algorithm
   * 
   * @param newAlg 
   */
  void setAlgorithm(innerOptAlgorithms newAlg)
  { alg = newAlg; }

  /** 
   * Limits of the hypercube. 
   * Currently, it assumes that all dimensions have the same limits.
   * 
   * @param down 
   * @param up 
   */
  void setLimits(double down, double up)
  {
    mDown = down;   mUp = up;
  }

  /** 
   * Compute the inner optimization algorithm
   * 
   * @param Xnext input: initial guess, output: result
   * 
   * @return error_code
   */
  int innerOptimize(vectord &Xnext);
  int innerOptimize(double* x, int n, void* objPointer);	


  /** 
   * Virtual function to be overriden by the actual function to be evaluated
   * 
   * @param query input point
   * 
   * @return function value at query point
   */
  virtual double innerEvaluate(const vectord& query) 
  {return 0.0;};


  /** 
   * Virtual function to be overriden by the actual function to be evaluated
   * 
   * @param query input point
   * @param grad output gradient at query point
   * 
   * @return function value at query point
   */
  virtual double innerEvaluate(const vectord& query, 
			       vectord& grad) 
  {return 0.0;};


protected:

  innerOptAlgorithms alg;
  double mDown, mUp;
};

#endif
