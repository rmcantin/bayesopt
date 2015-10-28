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

function x_out = langermann(x_in)
% Bounds [3,5]
    a = [3,5,2,1,7];
    b = [5,2,1,4,9];
    c = [1,2,5,2,3];
    
    x = x_in(1); 
    y = x_in(2);
  
    t1 = exp(-(x-a).^2/pi -(y-b).^2/pi);
    t2 = cos(pi*(x-a).^2 + pi*(y-b).^2);
    
    x_out = sum(c.*t1.*t2);