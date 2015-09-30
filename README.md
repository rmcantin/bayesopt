BayesOpt: A Bayesian optimization library
=========================================

BayesOpt is an efficient implementation of the Bayesian optimization
methodology for nonlinear-optimization, experimental design,
stochastic bandits and hyperparameter tunning. In the literature it is
also called Sequential Kriging Optimization (SKO), Sequential
Model-Based Optimization (SMBO) or Efficient Global Optimization
(EGO).

The online HTML version of these docs:
<http://rmcantin.bitbucket.org/html/>

Bayesian optimization uses a distribution over functions to build a
model of the unknown function for we are looking the extrema, and then
apply some active learning strategy to select the query points that
provides most potential interest or improvement. Thus, it is a
sampling efficient method for nonlinear optimization, design of
experiments or bandits-like problems.


Getting and installing BayesOpt
-------------------------------

The library can be download from any of this sources:

- Download: <https://bitbucket.org/rmcantin/bayesopt>
- Mirror: <https://github.com/rmcantin/bayesopt>
- Mirror: <http://mloss.org/software/view/453/>

You can also get the *cutting-edge* version from the repositories:

    >> hg clone https://bitbucket.org/rmcantin/bayesopt

or the git mirror:

    >> git clone https://github.com/rmcantin/bayesopt



The install guide and documentation for Windows, Linux and MacOS:
- [Online install guide](http://rmcantin.bitbucket.org/html/install.html) 
  or [Local install guide](@ref install)


Getting involved
----------------

The best place to ask questions and discuss about BayesOpt is the
[bayesopt-discussion mailing
list](https://groups.google.com/forum/#!forum/bayesopt-discussion). Alternatively,
you may directly contact Ruben Martinez-Cantin <rmcantin@unizar.es>.

Please file bug reports or suggestions at:
https://bitbucket.org/rmcantin/bayesopt/issues


Using BayesOpt for academic or commercial purposes
--------------------------------------------------

BayesOpt is licensed under the GPL and it is free to use. However,
please consider these recomentations when using BayesOpt:

- If you use BayesOpt in a work that leads to a **scientific
publication**, we would appreciate it if you would kindly cite BayesOpt
in your manuscript. Cite BayesOpt as:

> Ruben Martinez-Cantin, **BayesOpt: A Bayesian Optimization
> Library for Nonlinear Optimization, Experimental Design and
> Bandits**. Journal of Machine Learning Research, 15(Nov):3735--3739, 2014.

The paper can be found at http://jmlr.org/papers/v15/martinezcantin14a.html

- In addition, if you **use a specific algorithm** (REMBO, GP-Hedge,
etc.), please also cite the corresponding work. The reference for each
specific algorithm can be found in the documentation.

- If you are using the library for **research or academic purposes** or to
build free software, please send an email to <rmcantin@unizar.es> with
a brief description or link to your interest for this code (one or two
lines). There will be a section with links to research/papers/software
that use BayesOpt.

- Commercial applications may also acquire a **commercial license**. Please
contact <rmcantin@unizar.es> for details.


----------------------------------------------------------------------

Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>

BayesOpt is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BayesOpt is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with BayesOpt. If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/78838f26cf7d3fa7f8bda3a05eb63017 "githalytics.com")](http://githalytics.com/rmcantin/bayesopt)