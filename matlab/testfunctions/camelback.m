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

function y = camelback(x)
% -2 < x(1) < 2, -1 < x(2) < 1
% 2 Global minima: f(x) = -1.0316 (0.0898, -0.7126), (-0.0898, 0.7126)
    tmp1 = (4 - 2.1 * x(1)^2 + (x(1)^4)/3) * x(1)^2;
    tmp2 = x(1)*x(2);
    tmp3 = (-4 + 4 * x(2)^2) * x(2)^2;
    y = tmp1 + tmp2 + tmp3;
end
