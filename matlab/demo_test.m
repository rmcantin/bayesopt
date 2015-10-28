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
clear all, close all
addpath('testfunctions')

params.n_iterations = 190;
params.n_init_samples = 10;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1];
params.kernel_hp_std = [10];
params.verbose_level = 1;
params.log_filename = 'matbopt.log';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Continuous optimization');
fun = 'branin'; n = 2;
lb = zeros(n,1);
ub = ones(n,1);

tic;
bayesoptcont(fun,n,params,lb,ub)
toc;
disp('Press INTRO');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Discrete optimization');
% We still use branin
% The set of points must be numDimension x numPoints.
np = 100;
xset = repmat((ub-lb),1,np) .* rand(n,np) - repmat(lb,1,np);

tic;
bayesoptdisc(fun, xset, params)
toc;
yset = zeros(np,1);
for i=1:np
    yset(i) = feval(fun,xset(:,i));
end;
[y_min,id] = min(yset);
disp('Actual optimal');
disp(xset(:,id));
disp(y_min);
disp('Press INTRO');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Continuous optimization');
fun = 'hartmann'; n = 6;
lb = zeros(n,1);
ub = ones(n,1);

tic;
bayesoptcont(fun,n,params,lb,ub)
toc;
