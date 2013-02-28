% 
% -------------------------------------------------------------------------
%    This file is part of BayesOpt, an efficient C++ library for 
%    Bayesian optimization.
%
%    Copyright (C) 2011-2012 Ruben Martinez-Cantin <rmcantin@unizar.es>
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
clear all, close all
addpath('testfunctions')

params.n_iterations = 100;
params.n_init_iterations = 50;
params.c_name = 'EI';
params.s_name = 'GAUSSIAN_PROCESS';
params.noise = 0.005;
params.kernel = 'MATERN_ISO3';
params.theta = [0.5];
params.verbose_level = 5;
params.log_filename = 'matbopt.log';
%params.k_s_name = 'kMaternISO1';

n = 5;

lb = ones(n,1)*pi/2;
ub = ones(n,1)*pi;

tic;
bayesopt('michalewicz',n,params,lb,ub)
toc;

% The set of points must be nDim x nPoints.
xset = rand(n,100);

tic;
bayesoptdisc('quadratic',xset, params);
toc;