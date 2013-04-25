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

#include <valarray>
#include "bayesoptcont.hpp"

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
    return sqr(x-0.5) + sin(20*x)*0.2;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; }

  void printOptimal()
  {
    std::cout << "Optimal:" << 0.5 << std::endl;
  }

};

#include "unistd.h"
using namespace std;
using namespace bayesopt;
#include "matplotpp.h"

int is_run=1;
size_t state_ii = 0;
BayesOptBase* GLOBAL_MODEL;
vector<double> lx,ly;

class MP :public MatPlot{ 
void DISPLAY(){
    // Create test data 
    int n=400;
    vector<double> x,y,z,su,sl,c;
    x=linspace(0,1,n);
    y = x; z = x; su = x; sl = x; c= x;
    vectord q(1);
    for(int i=0;i<n;++i)
      {
	q(0) = x[i];
	ProbabilityDistribution* pd = GLOBAL_MODEL->getSurrogateModel()->prediction(q);
	y[i] = pd->getMean();
	su[i] = y[i] + 3*pd->getStd();
	sl[i] = y[i] - 3*pd->getStd();
	c[i] = -GLOBAL_MODEL->evaluateCriteria(q);
	z[i] = GLOBAL_MODEL->evaluateSample(q);
      }
    //plot
    subplot(2,1,1);
    title("press r to run and stop");
    plot(x,y);
    plot(lx,ly);set("k");set("*");
    plot(x,su);set("g");
    plot(x,sl);set("g");
    plot(x,z);set("r");
   
    subplot(2,1,2);
    plot(x,c);

}
}mp;

void display(){mp.display(); }
void reshape(int w,int h){ mp.reshape(w,h); }
void idle( void )
{
  if(is_run)
    {
      if (state_ii < 300)
	++state_ii;
      GLOBAL_MODEL->stepOptimization(state_ii); 
      vectord last(1);
      double res = GLOBAL_MODEL->getSurrogateModel()->getLastSample(last);
      ly.push_back(res);
      lx.push_back(last(0));
    
    }
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y ){ mp.mouse(button,state,x,y); }
void motion(int x, int y ){mp.motion(x,y); }
void passive(int x, int y ){mp.passivemotion(x,y); }
void keyboard(unsigned char key, int x, int y){
    mp.keyboard(key, x, y); 
    if(key=='r'){ if(is_run==0){is_run=1;}else{is_run=0;}}
}

int main(int nargs, char *args[])
{
  size_t dim = 1;
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 10;
  parameters.n_iterations = 100;
  parameters.surr_name = S_GAUSSIAN_PROCESS;
  // parameters.kernel.hp_mean[0] = 0.5;
  // parameters.kernel.hp_std[0] = 100.0;
  // parameters.kernel.n_hp = 1;
  // parameters.mean.name = "mZero";
  //  parameters.crit_name = "cHedge(cEI,cLCB,cExpReturn,cOptimisticSampling)";
  // parameters.epsilon = 0.0;

  state_ii = 0;

  ExampleOneD* opt = new ExampleOneD(dim,parameters);
  vectord result(dim);
  GLOBAL_MODEL = opt;
  opt->initializeOptimization();
  
  glutInit(&nargs, args);
  glutCreateWindow(100,100,500,400);
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    

  delete opt;
  return 0;
}
