
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
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

#include "randgen.hpp"
#include "lhs.hpp"
#include "gridsampling.hpp"
#include "log.hpp"
#include "bayesopt.hpp"

namespace bayesopt
{
  
  // DiscreteModel::DiscreteModel(const vecOfvec &validSet):
  //   BayesOptBase(), mInputSet(validSet)
  // {} // Constructor


  DiscreteModel::DiscreteModel( const vecOfvec &validSet, 
				bopt_params parameters):
    BayesOptBase(validSet[0].size(),parameters), mInputSet(validSet)
  {    
    mDims = mInputSet[0].size();    
  } // Constructor

  DiscreteModel::DiscreteModel(const vectori &categories, 
			       bopt_params parameters):
   BayesOptBase(categories.size(),parameters)
  {    
    mDims = categories.size();    
    utils::buildGrid(categories,mInputSet);
  }


  DiscreteModel::~DiscreteModel()
  {} // Default destructor


  vectord DiscreteModel::getFinalResult()
  {
    return getPointAtMinimum();
  }

  vectord DiscreteModel::samplePoint()
  {   
    randInt sample(mEngine, intUniformDist(0,mInputSet.size()-1));
    return mInputSet[sample()];
  };

  double DiscreteModel::evaluateSampleInternal( const vectord &query )
  { return evaluateSample(query); }; 


  
  void DiscreteModel::plotStepData(size_t iteration, const vectord& xNext,
				   double yNext)
  {
    if(mParameters.verbose_level >0)
      { 
	FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
			  << mParameters.n_iterations << " | Total samples: " 
			  << iteration+1+mParameters.n_init_samples ;
	FILE_LOG(logINFO) << "Trying point at: " << xNext ;
	FILE_LOG(logINFO) << "Current outcome: " << yNext ;
	FILE_LOG(logINFO) << "Best found at: " << getPointAtMinimum() ; 
	FILE_LOG(logINFO) << "Best outcome: " <<  getValueAtMinimum() ;    
      }
  }


  void DiscreteModel::sampleInitialPoints(matrixd& xPoints, vectord& yPoints)
  {

    vecOfvec perms = mInputSet;
    
    // By using random permutations, we guarantee that 
    // the same point is not selected twice
    utils::randomPerms(perms,mEngine);
    
    // vectord xPoint(mInputSet[0].size());
    for(size_t i = 0; i < yPoints.size(); i++)
      {
	const vectord xP = perms[i];
	row(xPoints,i) = xP;
	yPoints(i) = evaluateSample(xP);
      }
  }
  

  void DiscreteModel::findOptimal(vectord &xOpt)
  {
    xOpt = *mInputSet.begin();
    double min = evaluateCriteria(xOpt);
    
    for(vecOfvecIterator it = mInputSet.begin();
	it != mInputSet.end(); ++it)
      {
	double current = evaluateCriteria(*it);
	if (current < min)
	  {
	    xOpt = *it;  
	    min = current;
	  }
      }
  }

}  // namespace bayesopt


