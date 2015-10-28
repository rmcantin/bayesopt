% BAYESOPTDISC Optimization (minimization) of discrete target function 
% using Bayesian optimization.
%
% Usage: [xmin, fmin] = bayesoptdisc(@function_handler, validset, params)
%        [xmin, fmin] = bayesoptdisc('function_name', validset, params)
%
%
% params is a struct which have the same fields as the C/C++ interface 
%   (see include/parameters.h)
%
% validset is the set of discrete points for discrete optimization,
%      stacked in a single matrix. Thus, it must be a d x n matrix.
%
% 
% -------------------------------------------------------------------------
%    This file is part of BayesOpt, an efficient C++ library for 
%    Bayesian optimization.
%
%    Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
%
%    BayesOpt is free software: you can redistribute it and/or modify it 
%    under the terms of the GNU Affero General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BayesOpt is distributed in the hope that it will be useful, but 
%    WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.
%
%    You should have received a copy of the GNU Affero General Public License
%    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------
%

