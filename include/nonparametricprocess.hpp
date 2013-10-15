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
#include "ublas_extra.hpp"
#include "kernel_functors.hpp"
#include "mean_functors.hpp"
#include "prob_distribution.hpp"
#include "inneroptimization.hpp"

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
  class BAYESOPT_API NonParametricProcess
  {
  public:
    NonParametricProcess(size_t dim, bopt_params parameters);
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
    int fitInitialSurrogate(bool learnTheta = true);
  
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


    // Getters and setters
    void setSamples(const matrixd &x, const vectord &y);
    void addSample(const vectord &x, double y);
    double getSample(size_t index, vectord &x);
    double getLastSample(vectord &x);
    inline vectord getPointAtMinimum() { return mGPXX[mMinIndex]; };
    inline double getValueAtMinimum() { return mGPY(mMinIndex); };
    inline size_t getNSamples() { return mGPY.size(); };
    
    virtual double getSignalVariance() = 0;
  

    /** Sets the kind of learning methodology for kernel hyperparameters */
    inline void setLearnType(learning_type l_type) { mLearnType = l_type; };


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

    /** Wrapper of setMean for the C structure */
    int setMean (mean_parameters mean, size_t dim);

    /** 
     * \brief Computes the score (eg:likelihood) of the kernel
     * parameters.
     * @param query set of parameters.
     * @return score
     */
    double evaluateKernelParams(const vectord& query);


  protected:
    /** 
     * \brief Computes the negative log likelihood of the data for all
     * the parameters.
     * @return value negative log likelihood
     */
    virtual double negativeTotalLogLikelihood() = 0;


    /** 
     * \brief Computes the negative log likelihood of the data for the
     * kernel hyperparameters.
     * @return value negative log likelihood
     */
    virtual double negativeLogLikelihood() = 0;

    /** 
     * \brief Precompute some values of the prediction that do not
     * depends on the query.
     * @return error code
     */
    virtual int precomputePrediction() = 0;


    /** 
     * Computes the Cholesky decomposition of the Correlation (Kernel
     * or Gram) matrix 
     * @return error code
     */
    int computeCholeskyCorrelation();

    /** 
     * Adds a new point to the Cholesky decomposition of the Correlation 
     * matrix.
     * @return error code
     */
    int addNewPointToCholesky(const vectord& correlation,
			      double selfcorrelation);


    int computeCorrMatrix(matrixd& corrMatrix);
    matrixd computeCorrMatrix();
    matrixd computeDerivativeCorrMatrix(int dth_index);
    vectord computeCrossCorrelation(const vectord &query);
    double computeSelfCorrelation(const vectord& query);

  private:
    /**
     * Computes the negative score of the data using cross validation.
     * @return negative score
     */
    double negativeCrossValidation();

    /** 
     * \brief Computes the negative log prior of the hyperparameters.
     * @return value negative log prior
     */
    double negativeLogPrior();

    inline void checkBoundsY( size_t i )
    {
      if ( mGPY(mMinIndex) > mGPY(i) )       mMinIndex = i;
      else if ( mGPY(mMaxIndex) < mGPY(i) )  mMaxIndex = i;
    };


  protected:
    vecOfvec mGPXX;                                              ///< Data inputs
    vectord mGPY;                                                ///< Data values
    
    matrixd mFeatM;           ///< Value of the mean features at the input points
    vectord mMu;                 ///< Mean of the parameters of the mean function
    vectord mS_Mu;    ///< Variance of the params of the mean function W=mS_Mu*I

    boost::scoped_ptr<ParametricFunction> mMean;    ///< Pointer to mean function

    matrixd mL;             ///< Cholesky decomposition of the Correlation matrix
    size_t dim_;
    learning_type mLearnType;
    KernelModel mKernel;

  private:
    const double mRegularizer;   ///< Std of the obs. model (also used as nugget)
    size_t mMinIndex, mMaxIndex;	
    MeanFactory mPFactory;

    //TODO: might be unnecesary
    vectord mMeanV;                           ///< Mean value at the input points

    InnerOptimization* kOptimizer;
  };

  /**@}*/
  
} //namespace bayesopt

#endif
