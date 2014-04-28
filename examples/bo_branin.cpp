/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2013 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include "testfunctions.hpp"

int main(int nargs, char *args[])
{
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = 100;
  par.use_random_seed = 0;
  par.verbose_level = 1;
  
  ExampleBranin branin(2,par);
  vectord result(2);

  branin.optimize(result);
  std::cout << "Result: " << result << "->" 
	    << branin.evaluateSample(result) << std::endl;
  branin.printOptimal();

  return 0;
}
