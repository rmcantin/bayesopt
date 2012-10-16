% BAYESOPT Optimization (minimization) of target function using bayesian
% optimization.
%
% Usage: [xmin, fmin] = bayesopt(@function_handler, nDimensions, params)
%        [xmin, fmin] = bayesopt(@function_handler, nDimensions, params,
%                                lowerBound, upperBound) 
%
%        [xmin, fmin] = bayesopt('function_name', nDimensions, params)
%        [xmin, fmin] = bayesopt('function_name', nDimensions, params,
%                                lowerBound, upperBound) 
%
%
% nDimensions is the number of dimensions of the query vector.
%
% Params is a struct which may have the fields 
%      n_iterations           Maximum BayesOpt evaluations (budget) 
%      n_inner_iterations     Maximum inner optimizer evaluations 
%      n_init_samples         Number of samples before optimization 
%      verbose_level          Verbose level 
%      theta                  Vector with kernel hyperparameters 
%      alpha, beta, delta     Inv-Gamma-Normal hyperparameters 
%      noise                  Observation noise 
%      s_name            Name of the surrogate function 
%      k_name            Name of the kernel function 
%      c_name            Name of the criterion 
%      m_name            Name of the mean function 
%
% lowerBound and upperBound should be a nDim x 1 or 1 x nDim vectors with
%      the lower and upper bound for each component. (optional, default 0-1)
%
%
%-----------------------------------------------------------------------------
%   This file is part of BayesOptimization, an efficient C++ library for 
%   Bayesian optimization.
%
%   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
% 
%   BayesOptimization is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   BayesOptimization is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

