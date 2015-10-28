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

% BayesOpt parameters
params.n_iterations = 200;
params.n_init_samples = 10;
params.crit_name = 'cEI';
params.surr_name = 'sStudentTProcessNIG';
params.noise = 1e-6;
params.kernel_name = 'kMaternARD5';
params.kernel_hp_mean = [1 1];
params.kernel_hp_std = [10 10];
params.verbose_level = 5;
params.log_filename = 'matbopt.log';

fun = 'braninhighdim';    % the function has an effective 2D
n = 2;                    % number of low dims (effective)
nh = 1000;                % number of actual dims

global MATRIX_A            % random proyection matrix
global TRUE_I               % indexes of relevant components

TRUE_I = [150,237];

% Bounds of proyected (low-dim) space
lb = ones(n,1)*-sqrt(n);
ub = ones(n,1)*sqrt(n);

% rembo iterations results
nrembo = 10;                  % number of rembo iterations
values = zeros(nrembo,1);  
points = zeros(nrembo,n);    

for i=1:nrembo
    disp('Running optimization');
    MATRIX_A = randn(nh,n);
    
    tic;
    [result, val, error] = bayesoptcont(fun,n,params,lb,ub);
    toc;
    
    % We project the result back to the real low-dim space
    hd_res = MATRIX_A*result';  
    points(i,:) = hd_res(TRUE_I)';
    values(i) = val;
    
    fprintf('Opt point: %f, %f ||', hd_res(TRUE_I)); 
    fprintf('Opt value: %f || Error: %d\n',val, error);
end;

[foo,id] = min(values);
fprintf('Final point: %f, %f ||', points(id,:));  
fprintf('Final value: %f\n',values(id));