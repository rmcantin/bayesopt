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

#include "mixture_distribution.hpp"

MixtureDistribution::MixtureDistribution(size_t n):
  ProbabilityDistribution()
{
  
};

MixtureDistribution::~MixtureDistribution(){};

double MixtureDistribution::pdf(double x)
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->pdf(x);
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::negativeExpectedImprovement(double min,size_t g)
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->negativeExpectedImprovement(min,g);
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::lowerConfidenceBound(double beta)
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->lowerConfidenceBound(beta);
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::negativeProbabilityOfImprovement(double min,
							     double epsilon)
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->negativeProbabilityOfImprovement(min,epsilon);
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::sample_query(randEngine& eng)
{
  size_t n = mW.size();
  vectord cumw(n);
  
  std::partial_sum(mW.begin(), mW.end(), cumw.begin(), 
		   std::plus<double>());
  
  randFloat sampleUniform(eng, realUniformDist(0,1));
  double u = sampleUniform();

  for (size_t i=0; i < cumw.size(); ++i)
    {
      if (u < cumw(i))
	return mPD[i]->sample_query(eng);
    }
  return mPD[0]->sample_query(eng); //just in case...
};

double MixtureDistribution::getMean()
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->getMean();
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::getStd()
{
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      res(i) = mPD[i]->getStd();
    }
  return inner_prod(mW,res);
};

double MixtureDistribution::getGaussianStd()
{
  double totalMean = getMean();
  size_t n = mW.size();
  vectord res(n);
  for(size_t i=0;i<n;++i)
    {
      double sigma = mPD[i]->getStd();
      double mm = mPD[i]->getMean() - totalMean; 
      res(i) = mm*mm + sigma; 
    }
  return inner_prod(mW,res) ; 
}
