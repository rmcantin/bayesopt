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

params.n_iterations = 300;
params.n_init_iterations = 50;
params.crit_name = 'cEI';
params.surr_name = 'sGaussianProcessNormal';
params.noise = 0.005;
params.kernel_name = 'kMaternISO3';
params.kernel_hp_mean = [0.5];
params.kernel_hp_std = [10];
params.verbose_level = 0;
params.log_filename = 'matbopt.log';

n = 2;          % number of low dims (effective)
nh = 1000;      % number of actual dims
nreembo = 5;    % number of reembo iterations


global MATRIX_A
global truei

truei = [150,237];

lb = ones(n,1)*-sqrt(n);
ub = ones(n,1)*sqrt(n);
fun = 'braninhighdim';    % the function has an effective 2D
values = zeros(nreembo,1);
points = zeros(nreembo,n);

for i=1:nreembo
    disp('Continuous optimization');
    MATRIX_A = randn(nh,n);
    tic;
    result = bayesopt(fun,n,params,lb,ub);
    toc;

    values(i) = braninhighdim(result);
    hd_res = MATRIX_A*result';
    points(i,:) = hd_res(truei)';
    disp(hd_res(truei)); disp(values(i));
end;

[foo,id] = min(values);
disp('Final result');
disp(points(id,:));