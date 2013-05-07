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
%        [xmin, fmin] = bayesoptdisc(@function_handler, validset, params)
%        [xmin, fmin] = bayesoptdisc('function_name', validset, params)
%
%
% nDimensions is the number of dimensions of the query vector.
%
% Params is a struct which have the same fields as the C/C++ interface 
%   (see include/parameters.h)
%
% lowerBound and upperBound should be a nDim x 1 or 1 x nDim vectors with
%      the lower and upper bound for each component. (optional, default 0-1)
%
% validset is the set of discrete points for discrete optimization
%
% 
% -------------------------------------------------------------------------
%    This file is part of BayesOpt, an efficient C++ library for 
%    Bayesian optimization.
%
%    Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
%
%    BayesOpt is free software: you can redistribute it and/or modify it 
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BayesOpt is distributed in the hope that it will be useful, but 
%    WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------
%

