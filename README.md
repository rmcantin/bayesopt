BayesOpt: A Bayesian optimization toolbox            {#mainpage}
=========================================

BayesOpt is an free, efficient, implementation of the Bayesian
optimization methodology for nonlinear-optimization, experimental
design and stochastic bandits. In the literature it is also called
Sequential Kriging Optimization (SKO) or Efficient Global
Optimization (EGO). 

The online HTML version of these docs:
<http://rmcantin.bitbucket.org/html/>

Bayesian optimization uses a distribution over functions to build a
metamodel of the unknown function for we are looking the extrema,
and then apply some active learning strategy to select the query
points that provides most potential interest for the seek. For that
reason, it has been traditionally intended for optimization of
expensive function. However, the efficiency of the library make it
also interesting for many types of functions.

It is intended to be both fast and clear for development and
research. At the same time, it does everything the "right way". For
example:

- latin hypercube sampling is used for the preliminary design step,
- extensive use of Cholesky decomposition and related tecniques to 
  improve numeral stability and reduce computational cost,
- kernels, criteria and parametric functions can be combined to 
  produce more advanced functions,
- etc.

The documentation can be found in:
- [Install guide](http://rmcantin.bitbucket.org/html/install.html) | \ref install
- [Reference manual](http://rmcantin.bitbucket.org/html/reference.html)  | \ref reference
- [Bayesian optimization](http://rmcantin.bitbucket.org/html/bopttheory.html) | \ref bopttheory
- [Models and functions](http://rmcantin.bitbucket.org/html/modelopt.html) | \ref modelopt

You can also find more details at the [proyect wiki]
(http://bitbucket.org/rmcantin/bayesopt/wiki/Home) or subscribe to the
[bayesopt-discussion mailing
list](https://groups.google.com/forum/#!forum/bayesopt-discussion).

**Important:** This code is free to use. However, if you are using the
library, specially if it is for research or academic purposes, please
send me an email at <rmcantin@unizar.es> with your name, institution
and a brief description of your interest for this code (one or two
lines).

If you use BayesOpt in a work that leads to a scientific publication,
we would appreciate it if you would kindly cite BayesOpt in your
manuscript. Cite BayesOpt as something like:

----------------------------------------------------------------------

Ruben Martinez-Cantin, **BayesOpt: a toolbox for
nonlinear-optimization, experimental design and stochastic bandits**,
<http://bitbucket.org/rmcantin/bayesopt>

----------------------------------------------------------------------

Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>

BayesOpt is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BayesOpti is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with BayesOpt. If not, see <http:www.gnu.org/licenses/>.

----------------------------------------------------------------------