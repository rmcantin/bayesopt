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
  direct, lbfgs, combined
};


class InnerOptimization
{
public:
  InnerOptimization()
  { 
    alg = direct;    mDown = 0.;    mUp = 1.;
  };

  virtual ~InnerOptimization(){};

  virtual double innerEvaluate(const vectord& query, 
			       vectord& grad)
  {return 0.0;}

  void setAlgorithm(innerOptAlgorithms newAlg)
  { alg = newAlg; }

  void setLimits(double down, double up)
  {
    mDown = down;   mUp = up;
  }


protected:

  int innerOptimize(vectord &Xnext);
  int innerOptimize( double* x, int n, void* objPointer);	

  innerOptAlgorithms alg;
  double mDown, mUp;
};

#endif
