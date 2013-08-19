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
clear all, close all
addpath('testfunctions')

params.n_iterations = 100;
params.n_init_samples = 5;
params.crit_name = 'cEI';
params.surr_name = 'sGaussianProcessNormal';
params.noise = 0.005;
params.kernel_name = 'kMaternISO3';
params.kernel_hp_mean = [0.5];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';

% n = 5;
% lb = ones(n,1)*pi/2;
% ub = ones(n,1)*pi;
% fun = 'michalewicz';

n = 2;
lb = zeros(n,1);
ub = ones(n,1);
fun = 'branin';

disp('Continuous optimization');
tic;
bayesoptcont(fun,n,params,lb,ub)
toc;

disp('Discrete optimization');
% The set of points must be numDimension x numPoints.
np = 100;
xset = repmat((ub-lb),1,np) .* rand(n,np) - repmat(lb,1,np);

tic;
bayesoptdisc(fun,xset, params);
toc;