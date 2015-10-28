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

function y = hartmann(x)
%Bounds [0,1]^2
% Min = (0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573),
%        => -3.32236
    a = [10.0,   3.0, 17.0,   3.5,  1.7,  8.0;
          0.05, 10.0, 17.0,   0.1,  8.0, 14.0;
          3.0,   3.5,  1.7,  10.0, 17.0,  8.0;
         17.0,   8.0,  0.05, 10.0,  0.1, 14.0];
    c = [1.0, 1.2, 3.0, 3.2];
    p = [0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886;
         0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991;
         0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650;
         0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381];
    
    y = 0;
    for i=1:4
        sum = 0;
        for j=1:6
            sum = sum - a(i,j)*(x(j)-p(i,j))^2;
        end;
        y = y - c(i)*exp(sum);
    end;
