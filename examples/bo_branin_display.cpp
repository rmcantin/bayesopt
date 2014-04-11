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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include "bayesopt.hpp"
#include "displaygp.hpp"

class ExampleBranin: public bayesopt::ContinuousModel
{
public:
  ExampleBranin(size_t dim,bopt_params par):
    ContinuousModel(dim,par) {}

  double evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }

    double x = xin(0) * 15 - 5;
    double y = xin(1) * 15;

    return sqr(y-(5.1/(4*sqr(M_PI)))*sqr(x)+5*x/M_PI-6)+10*(1-1/(8*M_PI))*cos(x)+10;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(2);  
    sv(0) = 0.1239; sv(1) = 0.8183;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.5428; sv(1) = 0.1517;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.9617; sv(1) = 0.1650;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

};

// Unfortunately OpenGL functions require no parameters, so the object has to be global.
bayesopt::utils::DisplayProblem2D GLOBAL_MATPLOT;

void display( void ){ GLOBAL_MATPLOT.display(); }
void reshape( int w,int h ){ GLOBAL_MATPLOT.reshape(w,h); }
void idle( void ) { glutPostRedisplay(); } 

void mouse(int button, int state, int x, int y ){ GLOBAL_MATPLOT.mouse(button,state,x,y); }
void motion(int x, int y ){ GLOBAL_MATPLOT.motion(x,y); }
void passive(int x, int y ){ GLOBAL_MATPLOT.passivemotion(x,y); }

void keyboard(unsigned char key, int x, int y)
{
    GLOBAL_MATPLOT.keyboard(key, x, y); 
    if(key=='r')   //Toogle run/stop
      { 
	GLOBAL_MATPLOT.toogleRUN();
      }
    if(key=='s')   //Activate one step
      { 
	GLOBAL_MATPLOT.setSTEP();
      }
}


int main(int nargs, char *args[])
{
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = 100;
  par.n_inner_iterations = 1000;
  par.n_init_samples = 2;
  par.n_iter_relearn = 1;
  par.use_random_seed = 0;
  
  par.l_all = 0;
  par.l_type = L_MCMC;
  par.sc_type = SC_MAP;
  par.verbose_level = 1;

  par.kernel.name = "kMaternARD5";
  par.kernel.hp_mean[0] = 1.0;
  par.kernel.hp_std[0] = 10.0;
  par.kernel.hp_mean[1] = 1.0;
  par.kernel.hp_std[1] = 10.0;
  par.kernel.n_hp = 2;

  par.surr_name = "sStudentTProcessNIG";
  par.sigma_s = 10;
  par.noise = 0.01;

  // par.mean.name = new char[128];
  // strcpy(par.mean.name,"mConst");
  // par.mean.coef_mean[0] = 50;
  // par.mean.coef_std[0] = 10;
  // par.mean.n_coef = 1;
  
  // par.crit_name = "cHedge(cLCB,cEI,cPOI)";
  // par.n_iter_relearn = 20;
  // double cParams[] = {5.0, 1.0, 0.01};
  // std::copy(cParams, cParams+3, par.crit_params);
  // par.n_crit_params = 3;
  
  boost::scoped_ptr<ExampleBranin> branin(new ExampleBranin(2,par));
  GLOBAL_MATPLOT.init(branin.get(),2);

  vectord sv(2);  
  sv(0) = 0.1239; sv(1) = 0.8183;
  GLOBAL_MATPLOT.setSolution(sv);
  
  sv(0) = 0.5428; sv(1) = 0.1517;
  GLOBAL_MATPLOT.setSolution(sv);

  sv(0) = 0.9617; sv(1) = 0.1650;
  GLOBAL_MATPLOT.setSolution(sv);

  glutInit(&nargs, args);
  glutCreateWindow(50,50,800,650);
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    


  // vectord result(2);

  // branin.optimize(result);
  // std::cout << "Result: " << result << "->" 
  // 	    << branin.evaluateSample(result) << std::endl;
  // branin.printOptimal();

  return 0;
}
