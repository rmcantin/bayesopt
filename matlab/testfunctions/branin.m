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

function y = branin(x)
%Bounds [0,1]^2
% Min = 0.1239 0.8183
% Min = 0.5428 0.1517  => 0.3979
% Min = 0.9617 0.1650
    
a = x(1) * 15 - 5;
b = x(2) * 15;

y = (b-(5.1/(4*pi^2))*a^2+5*a/pi-6)^2+10*(1-1/(8*pi))*cos(a)+10;
