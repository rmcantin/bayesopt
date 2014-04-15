/**  \file displaygp.hpp 
     \brief Plots the evolution (nonparametric process and criteria)
     of a 1D problem. */
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

#ifndef _DISPLAYGP_HPP_
#define _DISPLAYGP_HPP_

//#include <stdexcept>
#include "bayesopt.hpp"
#include "matplotpp.h"  

namespace bayesopt
{
  namespace utils
  {      
    enum RunningStatus
      {
	RUN, STEP, STOP, NOT_READY
      };

    class DisplayProblem1D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      BayesOptBase* bopt_model;
      std::vector<double> lx,ly;

    public:
      DisplayProblem1D(): MatPlot()
      {
	status = NOT_READY;
      }

      void init(BayesOptBase* bopt, size_t dim)
      {
	if (dim != 1) 
	  { 
	    throw std::invalid_argument("Display only works for 1D problems"); 
	  }

	bopt_model = bopt;
	bopt->initializeOptimization();
	size_t n_points = bopt->getSurrogateModel()->getData()->getNSamples();
	for (size_t i = 0; i<n_points;++i)
	  {
	    const double res = bopt->getSurrogateModel()->getData()->getSampleY(i);
	    const vectord last = bopt->getSurrogateModel()->getData()->getSampleX(i);
	    ly.push_back(res);
	    lx.push_back(last(0));
	  }
	state_ii = 0;    
	status = STOP;
      };

      void setSTEP()
      {
	if (status != NOT_READY)
	  {
	    status = STEP;
	  }
      };

      void toogleRUN()
      {
	if (status != NOT_READY)
	  {
	    if(status != RUN)
	      {
		status = RUN;
	      }
	    else
	      {
		status = STOP;
	      }
	  }
      }

      void DISPLAY()
      {
	if (status != NOT_READY)
	  {
	    size_t nruns = bopt_model->getParameters()->n_iterations;
	    if ((status != STOP) && (state_ii < nruns))
	      {
		// We are moving. Next iteration
		++state_ii;
		bopt_model->stepOptimization(state_ii); 
		const double res = bopt_model->getSurrogateModel()->getData()->getLastSampleY();
		const vectord last = bopt_model->getSurrogateModel()->getData()->getLastSampleX();
		ly.push_back(res);
		lx.push_back(last(0));
	  
		if (status == STEP) { status = STOP; }
	      }

	    // We compute the prediction, true value and criteria at 1000 points
	    int n=1000;
	    std::vector<double> x,y,z,su,sl,c;
	    x = linspace(0,1,n);
	    y = x; z = x; su = x; sl = x; c = x;

	    // Query functions at the 1000 points
	    vectord q(1);
	    for(size_t i=0; i<n; ++i)
	      {
		q(0) = x[i];                                                 // Query
		ProbabilityDistribution* pd = bopt_model->getSurrogateModel()->prediction(q);
		y[i] = pd->getMean();                                //Expected value
		su[i] = y[i] + 2*pd->getStd();                       //Upper bound (95 %)
		sl[i] = y[i] - 2*pd->getStd();                       //Lower bound (95 %)
		c[i] = -bopt_model->getCriteria()->evaluate(q);      //Criteria value
		z[i] = bopt_model->evaluateSample(q);                //Target function true value
	      }
 
	    //GP subplot
	    subplot(2,1,1);
	    title("Press r to run and stop, s to run a step and q to quit.");
	    plot(x,y); set(3);                            // Expected value in default color (blue)
	    plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
	    plot(x,su);set("g"); set(2);                  // Uncertainty as green lines
	    plot(x,sl);set("g"); set(2);
	    plot(x,z);set("r"); set(3);                   // True function as red line
	    
	    //Criterion subplot
	    subplot(2,1,2);
	    plot(x,c); set(3);
	  }
      };
    };

    class DisplayProblem2D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      BayesOptBase* bopt_model;
      std::vector<double> lx,ly;
      std::vector<double> cx, cy;
      std::vector<double> solx, soly;

    public:
      DisplayProblem2D(): 
	MatPlot(), cx(1), cy(1)
      {
	status = NOT_READY;
      }

      void setSolution(vectord sol)
      {
	solx.push_back(sol(0));
	soly.push_back(sol(1));
      }

      void init(BayesOptBase* bopt, size_t dim)
      {
	if (dim != 2) 
	  { 
	    throw std::invalid_argument("This display only works for 2D problems"); 
	  }

	bopt_model = bopt;
	bopt->initializeOptimization();
	size_t n_points = bopt->getSurrogateModel()->getData()->getNSamples();
	for (size_t i = 0; i<n_points;++i)
	  {
	    const vectord last = bopt->getSurrogateModel()->getData()->getSampleX(i);
	    lx.push_back(last(0));
	    ly.push_back(last(1));
	  }
	state_ii = 0;    
	status = STOP;
      };

      void setSTEP()
      {
	if (status != NOT_READY)
	  {
	    status = STEP;
	  }
      };

      void toogleRUN()
      {
	if (status != NOT_READY)
	  {
	    if(status != RUN)
	      {
		status = RUN;
	      }
	    else
	      {
		status = STOP;
	      }
	  }
      }

      void DISPLAY()
      {
	if (status != NOT_READY)
	  {
	    size_t nruns = bopt_model->getParameters()->n_iterations;
	    title("Press r to run and stop, s to run a step and q to quit.");
	    plot(cx,cy);set("g");set("o");set(4);         // Data points as black star
	    plot(solx,soly);set("r"); set("o");set(4);    // Solutions as red points

	    if ((status != STOP) && (state_ii < nruns))
	      {
		// We are moving. Next iteration
		++state_ii;
		bopt_model->stepOptimization(state_ii); 
		const vectord last = bopt_model->getSurrogateModel()->getData()->getLastSampleX();
		//GP subplot
		cx[0] = last(0);
		cy[0] = last(1);

		if (!lx.empty())
		  {	
		    plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
		  }

		lx.push_back(last(0));
		ly.push_back(last(1));

		if (status == STEP) { status = STOP; }
	      }	    
	    else
	      {
		plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
	      }

	  }
      };
    };


  } //namespace utils

} //namespace bayesopt


#endif
