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
% MATLAB
if (ispc)
    if exist('../bin/Release/bayesopt.dll','file')
        disp('Compiling dynamic library');
        mex -DBAYESOPT_DLL -output bayesoptcont bayesoptmex.c ...
            -L..\lib\Release -L. -lbayesopt ...
            -I..\include -I..\wrappers
        mex -DBAYESOPT_DLL -output bayesoptdisc bayesoptdiscmex.c ...
            -L..\lib\Release -L. -lbayesopt ...
            -I..\include -I..\wrappers
        mex -DBAYESOPT_DLL -output bayesoptcat bayesoptcatmex.c ...
            -L..\lib\Release -L. -lbayesopt ...
            -I..\include -I..\wrappers
    else
        disp('Compiling static library');
        mex -output bayesoptcont bayesoptmex.c  ...
            -L../lib/Release -lbayesopt -lnlopt ...
            -I../include -I../wrappers
        
        mex -output bayesoptdisc bayesoptdiscmex.c  ...
            -L../lib/Release -lbayesopt -lnlopt ...
            -I../include -I../wrappers

        mex -output bayesoptcat bayesoptcatmex.c  ...
            -L../lib/Release -lbayesopt -lnlopt ...
            -I../include -I../wrappers
    end
else
    mex -output bayesoptcont bayesoptmex.c -L../lib -lbayesopt ...
        -lnlopt -I../include -I../wrappers -I../nlopt/api 

    mex -output bayesoptdisc bayesoptdiscmex.c -L../lib -lbayesopt ...
        -lnlopt -I../include -I../wrappers -I../nlopt/api 

    mex -output bayesoptcat bayesoptcatmex.c -L../lib -lbayesopt ...
        -lnlopt -I../include -I../wrappers -I../nlopt/api 
end