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

% You can also change ../lib for the correspoding install path
% Octave
mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
    --mex --output bayesoptcont.mex bayesoptmex.c 

mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
    --mex --output bayesoptdisc.mex bayesoptdiscmex.c

mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
    --mex --output bayesoptcat.mex bayesoptcatmex.c
