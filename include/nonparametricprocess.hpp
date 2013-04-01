/** \file nonparametricprocess.hpp 
    \brief Nonparametric process abstract module */
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


#ifndef __NONPARAMETRICPROCESS_HPP__
#define __NONPARAMETRICPROCESS_HPP__

#include <boost/scoped_ptr.hpp>
#include <boost/math/distributions/normal.hpp> 
#include "parameters.h"
#include "kernel_functors.hpp"
#include "mean_functors.hpp"
#include "specialtypes.hpp"
#include "inneroptimization.hpp"	
#include "prob_distribution.hpp"

#define USE_CHOL 1

namespace bayesopt
{

  /** \addtogroup NonParametricProcesses 
   *  \brief Set of nonparametric processes (Gaussian process, Student's
   *  t process, etc.) for surrogate modelling
   */
  /**@{*/

  /**
   * \brief Abstract class to implement non-parametric processes
   */
  class NonParametricProcess: public InnerOptimization
  {
  public:
    NonParametricProcess(size_t dim, double noise);
    virtual ~NonParametricProcess();

    /** 
     * \brief Factory model generator for surrogate models
     * @param parameters (process name, noise, priors, etc.)
     * @return pointer to the corresponding derivate class (surrogate model)
     */
    static NonParametricProcess* create(size_t dim, bopt_params parameters);

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    virtual ProbabilityDistribution* prediction(const vectord &query) = 0;
		 		 
    /** 
     * \brief Computes the initial surrogate model. 
     * It also estimates the kernel parameters by MAP or ML. This
     * function is hightly inefficient.  Use it only at the begining.
     * @return error code
     */
    int fitInitialSurrogate();
  
    /** 
     * \brief Sequential update of the surrogate model by adding a new point.
     *  Add new point efficiently using Matrix Decomposition Lemma 
     *  for the inversion of the correlation matrix or adding new 
     *  rows to the Cholesky decomposition.
     * @return error code
     */   
    int updateSurrogateModel( const vectord &Xnew,
			      double Ynew);

    /** 
     * \brief Full update of the surrogate model by adding a new point.
     *  It recomputes the kernel hyperparameters and full covariance matrix.
     * @return error code
     */   
    int fullUpdateSurrogateModel( const vectord &Xnew,
				  double Ynew);


    /// Getters and setters
    void setSamples(const matrixd &x, const vectord &y);
    void addSample(const vectord &x, double y);
    double getSample(size_t index, vectord &x);
    double getLastSample(vectord &x)
    {
      size_t last = mGPY.size()-1;
      x = mGPXX[last];
      return mGPY[last];
    };
    inline vectord getPointAtMinimum() { return mGPXX[mMinIndex]; };
    inline double getValueAtMinimum() { return mGPY(mMinIndex); };
  

    // /** 
    //  * \brief Select kernel (covariance function) for the surrogate process.
    //  * @param thetav kernel parameters
    //  * @param k_name kernel name
    //  * @return error_code
    //  */
    // int setKernel (const vectord &thetav, kernel_name k_name, size_t dim);

    // /** 
    //  * \brief Wrapper of setKernel for c arrays
    //  */
    // inline int setKernel (const double *theta, size_t n_theta, 
    // 			  kernel_name k_name, size_t dim)
    // {
    //   vectord th(n_theta);
    //   std::copy(theta, theta+n_theta, th.begin());
    //   int error = setKernel(th, k_name, dim);
    // };

    /** 
     * \brief Select kernel (covariance function) for the surrogate process.
     * @param thetav kernel parameters (mean)
     * @param stheta kernel parameters (std)
     * @param k_name kernel name
     * @return error_code
     */
    int setKernel (const vectord &thetav, const vectord &stheta, 
		   std::string k_name, size_t dim);

    /** 
     * \brief Wrapper of setKernel for c arrays
     */
    inline int setKernel (kernel_parameters kernel, size_t dim)
    {
      size_t n = kernel.n_theta;
      vectord th(n);
      vectord sth(n);
      std::copy(kernel.theta, kernel.theta+n, th.begin());
      std::copy(kernel.s_theta, kernel.s_theta+n, sth.begin());
      int error = setKernel(th, sth, kernel.name, dim);
    };

    /** Set prior (Gaussian) for kernel hyperparameters */
    inline int setKernelPrior (const vectord &theta, const vectord &s_theta)
    {
      size_t n_theta = theta.size();
      for (size_t i = 0; i<n_theta; ++i)
	{
	  boost::math::normal n(theta(i),s_theta(i));
	  priorKernel.push_back(n);
	}
      return 0;
    };



    // /** 
    //  * \brief Select the parametric part of the surrogate process.
    //  * 
    //  * @param muv mean function parameters
    //  * @param m_name mean function name
    //  * @return error_code
    //  */
    // int setMean (const vectord &muv, mean_name m_name, size_t dim);

    // /** 
    //  * \brief Wrapper of setMean for c arrays
    //  */
    // inline int setMean (const double *mu, size_t n_mu, 
    // 			mean_name m_name, size_t dim)
    // {
    //   vectord vmu(n_mu);
    //   std::copy(mu, mu+n_mu, vmu.begin());
    //   return setMean(vmu, m_name, dim);
    // };


    /** 
     * \brief Select the parametric part of the surrogate process.
     * 
     * @param muv mean function parameters
     * @param smu std function parameters
     * @param m_name mean function name
     * @param dim number of input dimensions
     * @return error_code
     */
    int setMean (const vectord &muv, const vectord &smu, 
		 std::string m_name, size_t dim);

    /** 
     * \brief Wrapper of setMean for c arrays
     */
    inline int setMean (mean_parameters mean, size_t dim)
    {
      size_t n_mu = mean.n_mu;
      vectord vmu(n_mu);
      vectord smu(n_mu);
      std::copy(mean.mu, mean.mu+n_mu, vmu.begin());
      std::copy(mean.s_mu, mean.s_mu+n_mu, smu.begin());
      return setMean(vmu, smu, mean.name, dim);
    };


  protected:

    /** 
     * \brief Computes the negative log likelihood of the data.
     * @return value negative log likelihood
     */
    virtual double negativeLogLikelihood() = 0;

    /** 
     * \brief Precompute some values of the prediction that do not depends on
     * the query
     * @return error code
     */
    virtual int precomputePrediction() = 0;

    double innerEvaluate(const vectord& query)
    { 
      mKernel->setHyperParameters(query);
      double posterior = negativeLogLikelihood();
      for(size_t i = 0; i<query.size();++i)
	{
	  if (priorKernel[i].standard_deviation() > 0)
	    {
	      posterior *= pdf(priorKernel[i],query(i));
	    }
	}
      return posterior;
    };


    /** 
     * Computes the inverse of the Correlation (Kernel or Gram) matrix  
     * @return error code
     */
    int computeInverseCorrelation();
    int addNewPointToInverse(const vectord& correlation,
			     double selfcorrelation);

    int computeCholeskyCorrelation();
    int addNewPointToCholesky(const vectord& correlation,
			      double selfcorrelation);


    int computeCorrMatrix(matrixd& corrMatrix);
    matrixd computeCorrMatrix();
    matrixd computeDerivativeCorrMatrix(int dth_index);
    vectord computeCrossCorrelation(const vectord &query);

    inline void checkBoundsY( size_t i )
    {
      if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
      else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
    };


  protected:
    const double mRegularizer;   ///< Std of the obs. model (also used as nugget)
    vecOfvec mGPXX;                                              ///< Data inputs
    vectord mGPY;                                                ///< Data values
    vectord mMeanV;                           ///< Mean value at the input points
    matrixd mFeatM;           ///< Value of the mean features at the input points

    std::vector<boost::math::normal> priorKernel;   ///< Prior of kernel function
    boost::scoped_ptr<Kernel> mKernel;            ///< Pointer to kernel function
    boost::scoped_ptr<ParametricFunction> mMean;    ///< Pointer to mean function

    // TODO: Choose one
    matrixd mL;             ///< Cholesky decomposition of the Correlation matrix
    covMatrix mInvR;                              ///< Inverse Correlation matrix
    size_t dim_;

  private:
    size_t mMinIndex, mMaxIndex;	
    KernelFactory mKFactory;
    MeanFactory mPFactory;
  };

  /**@}*/
  
} //namespace bayesopt

#endif
