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

#include "bayesopt.hpp"
#include "displaygp.hpp"

class ExampleOneD: public bayesopt::ContinuousModel
{
public:
  ExampleOneD(size_t dim, bopt_params par):
    ContinuousModel(dim,par) {}

  double evaluateSample(const vectord& xin)
  {
    if (xin.size() > 1)
      {
	std::cout << "WARNING: This only works for 1D inputs." << std::endl
		  << "WARNING: Using only first component." << std::endl;
      }

    double x = xin(0);
    return (x-0.3)*(x-0.3) + sin(20*x)*0.2;
  };

  bool checkReachability(const vectord &query) {return true;};

  void printOptimal()
  {
    std::cout << "Optimal:" << 0.5 << std::endl;
  }

};


// Unfortunately OpenGL functions require no parameters, so the object has to be global.
bayesopt::utils::DisplayProblem1D GLOBAL_MATPLOT;

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
  size_t dim = 1;
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 7;
  parameters.n_iter_relearn = 0;
  parameters.n_iterations = 300;
  parameters.verbose_level = 2;

  // Surrogate models
  //  parameters.surr_name = "sStudentTProcessNIG";
  parameters.surr_name = "sGaussianProcessNormal";

  // Criterion model
  parameters.crit_name = "cAopt";
  parameters.n_crit_params = 0;

  // parameters.crit_name = "cEI";
  // parameters.crit_params[0] = 1;
  // parameters.n_crit_params = 1;

  // parameters.crit_name = "cLCB";
  // parameters.crit_params[0] = 5;
  // parameters.n_crit_params = 1;

  // Kernel models
  // parameters.kernel.name = "kSum(kPoly3,kRQISO)";
  // double mean[128] = {1, 1, 1, 1};
  // double std[128] = {10, 10, 10, 10};
  // size_t nhp = 4;
  // memcpy(parameters.kernel.hp_mean, mean, nhp * sizeof(double));
  // memcpy(parameters.kernel.hp_std,std, nhp * sizeof(double));
  // parameters.kernel.n_hp = nhp;

  parameters.kernel.name = "kMaternISO3";
  parameters.kernel.hp_mean[0] = 1;
  parameters.kernel.hp_std[0] = 5;
  parameters.kernel.n_hp = 1;

  boost::scoped_ptr<ExampleOneD> opt(new ExampleOneD(dim,parameters));
  GLOBAL_MATPLOT.init(opt.get(),1);

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

  return 0;
}
