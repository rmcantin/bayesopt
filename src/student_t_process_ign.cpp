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

#include "student_t_process.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"
#include "student_t_distribution.hpp"

using boost::numeric::ublas::inplace_solve;
using boost::numeric::ublas::lower_tag;
using boost::numeric::ublas::lower;

  
StudentTProcessIGN::StudentTProcessIGN(double noise, double alpha, 
					double beta, double delta):
  NonParametricProcess(noise),
  mAlpha(alpha), mBeta (beta), mDelta2(delta)
{}  // Constructor



StudentTProcessIGN::~StudentTProcessIGN()
{} // Default destructor


double StudentTProcessIGN::negativeLogLikelihoodIGN()
{
  
  vectord h = mMean->getFeatureVector()
