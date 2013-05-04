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


#include "cholesky.hpp"
#include "trace_ublas.hpp"
#include "student_t_process_jef.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas; 

  StudentTProcessJeffreys::StudentTProcessJeffreys(size_t dim, 
						   bopt_params params):
    HierarchicalGaussianProcess(dim, params)
  {
    d_ = new StudentTDistribution();
  }  // Constructor



  StudentTProcessJeffreys::~StudentTProcessJeffreys()
  {
    delete d_;
  } // Default destructor




  double StudentTProcessJeffreys::negativeLogLikelihood()
  {
    /* In this case, they are equivalent */
    return negativeTotalLogLikelihood();
  }


  ProbabilityDistribution* 
  StudentTProcessJeffreys::prediction(const vectord &query )
  {
    double kq = (*mKernel)(query, query);;
    vectord kn = computeCrossCorrelation(query);
    vectord phi = mMean->getFeatures(query);
  
    vectord v(kn);
    inplace_solve(mL,v,ublas::lower_tag());

    vectord rq = phi - prod(v,mKF);

    vectord rho(rq);
    inplace_solve(mL2,rho,ublas::lower_tag());
    
    double yPred = inner_prod(phi,mWML) + inner_prod(v,mAlphaF);
    double sPred = sqrt( mSigma * (kq - inner_prod(v,v) 
				   + inner_prod(rho,rho)));

    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }

  int StudentTProcessJeffreys::precomputePrediction()
  {
    size_t n = mGPXX.size();
    size_t p = mMean->nFeatures();

    mKF = trans(mFeatM);
    inplace_solve(mL,mKF,ublas::lower_tag());

    matrixd FKF = prod(trans(mKF),mKF);
    mL2 = FKF;
    utils::cholesky_decompose(FKF,mL2);

    vectord Ky(mGPY);
    inplace_solve(mL,Ky,ublas::lower_tag());

    mWML = prod(Ky,mKF);
    utils::cholesky_solve(mL2,mWML,ublas::lower());

    mAlphaF = mGPY - prod(mWML,mFeatM);
    inplace_solve(mL,mAlphaF,ublas::lower_tag());
    mSigma = inner_prod(mAlphaF,mAlphaF)/(n-p);
    
    d_->setDof(n-p);  
    return 1;
  }

} //namespace bayesopt
